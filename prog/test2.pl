use Parallel::ForkManager 0.7.6;
use Data::Dumper;  # to display the data structures retrieved.
use strict;
 
my $pm = Parallel::ForkManager->new(20);  # using the system temp dir $L<File::Temp::tempdir()
 
# data structure retrieval and handling
my %retrieved_responses = ();  # for collecting responses
$pm -> run_on_finish (
  sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
 
    # see what the child sent us, if anything
    if (defined($data_structure_reference)) {  # test rather than assume child sent anything
      my $reftype = ref($data_structure_reference);
      print qq|ident "$ident" returned a "$reftype" reference.\n\n|;
      if (1) {  # simple on/off switch to display the contents
        print &Dumper($data_structure_reference) . qq|end of "$ident" sent structure\n\n|;
      }
 
      # we can also collect retrieved data structures for processing after all children have exited
      $retrieved_responses{$ident} = $data_structure_reference;
    } else {
      print qq|ident "$ident" did not send anything.\n\n|;
    }
  }
);
 
# generate a list of instructions
my @instructions = (  # a unique identifier and what the child process should send
  {'name' => '%ENV keys as a string', 'send' => 'keys'},
  {'name' => 'Send Nothing'},  # not instructing the child to send anything back to the parent
  {'name' => 'Childs %ENV', 'send' => 'all'},
  {'name' => 'Child chooses randomly', 'send' => 'random'},
  {'name' => 'Invalid send instructions', 'send' => 'Na Na Nana Na'},
  {'name' => 'ENV values in an array', 'send' => 'values'},
);
 
my $instruction = '';
foreach $instruction (@instructions) {
  $pm->start($instruction->{'name'}) and next;  # this time we are using an explicit, unique child process identifier
 
  # last step in child processing
  $pm->finish(0) unless $instruction->{'send'};  # no data structure is sent unless this child is told what to send.
 
  if ($instruction->{'send'} eq 'keys') {
    $pm->finish(0, \join(', ', keys %ENV));
 
  } elsif ($instruction->{'send'} eq 'values') {
    $pm->finish(0, [values %ENV]);  # kinda useless without knowing which keys they belong to...
 
  } elsif ($instruction->{'send'} eq 'all') {
    $pm->finish(0, \%ENV);  # remember, we are not "returning" anything, just copying the hash to disc
 
  # demonstrate clearly that the child determines what type of reference to send
  } elsif ($instruction->{'send'} eq 'random') {
    my $string = q|I'm just a string.|;
    my @array = qw(I am an array);
    my %hash = (type => 'associative array', synonym => 'hash', cool => 'very :)');
    my $return_choice = ('string', 'array', 'hash')[int(rand 3)];  # randomly choose return data type
    $pm->finish(0, \$string) if ($return_choice eq 'string');
    $pm->finish(0, \@array) if ($return_choice eq 'array');
    $pm->finish(0, \%hash) if ($return_choice eq 'hash');
 
  # as a responsible child, inform parent that their instruction was invalid
  } else {
    $pm->finish(0, \qq|Invalid instructions: "$instruction->{'send'}".|);  # ordinarily I wouldn't include invalid input in a response...
  }
}
$pm->wait_all_children;  # blocks until all forked processes have exited
 
# post fork processing of returned data structures
for (sort keys %retrieved_responses) {
  print qq|Post processing "$_"...\n|;
}
