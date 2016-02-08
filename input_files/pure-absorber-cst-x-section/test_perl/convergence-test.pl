#!/usr/bin/perl

# the strict package forces you to declare each variable you use beforehand
use strict;

# a variable in strict mode is declared using my
# the $ symbol means it is a single-valued variable
# the @ symbol means it is an array
# each declaration/instruction is closed with a ; sign
my $cell_init = $ARGV[0];
my $nb_run = $ARGV[1];
my @cells_list = (0) x $nb_run;
@cells_list[0] = $cell_init;
 for (my $i=1; $i<$nb_run; $i++)
  {
    @cells_list[$i] = 2.*@cells_list[$i-1]
  }

# now, we read in a variable that will be the filename of the template input file
my $file_in = $ARGV[2];

## start of the loop
for( my $i=0; $i< $nb_run; $i++)
{
  my $value = @cells_list[$i];
  print "This is the current parameter value: $value \n";

  # now we create a new string variable that will later be the filename of the new input deck
  # the . symbol is the concatenation operator between strings
  my $new_input_filename = $file_in."-nb-cells-".$value.".i";
  print " The new filename is $new_input_filename \n";

  # open the template file and store its filehandle (fh_in)
  open my $fh_in,  '<', $file_in.".i" or die "Can't open output $file_in !";
  # open the new file (it currently does not exist and is thus empty) and store its filehandle (fh_out)
  open my $fh_out, '>', $new_input_filename or die "Can't open output $new_input_filename !";

  while (<$fh_in>)
  {
    # this is for you to see on the console, we read line-by-line, while there is something
    # the line read is stored in a special PERL variable $_
    # now we actually print that line into the new file
    # BUT BEFORE THAT, we change the dummy characters for the real value
    # we use a regular expression (read the tutorials for more details_
    # s = substitute
    my $param_to_change = "nx =";
    s/$param_to_change/$param_to_change$value/;
    $param_to_change = "file_base =";
    my $file_base = "file_base = ".$file_in."-nb-cells-".$value."_out";
    s/$param_to_change/$file_base/;
    print $fh_out $_;
  }
  close $fh_in;
  close fh_out;

  # run the input file

  system('./../../../rhea-opt -i '.$new_input_filename);
}

print " I am done with this !!! \n";
exit 111;
