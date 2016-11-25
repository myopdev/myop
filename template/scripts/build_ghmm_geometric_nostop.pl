#!/usr/bin/perl

use strict;
use warnings;
use File::Copy;
use Parallel::ForkManager;

my $metaparfile = "cnf/meta.cnf";
my $configdir = "ghmm/cnf/";
my $ncpu = 1;


my %metapar;
# read metaparameters file
open (META, "<$metaparfile") or die "Cant open $metaparfile: $!\n";
foreach my $line (<META>)
  {
    chomp($line);
    my @fields = split(/\s*=\s*/, $line);
    # remove spaces;
    $fields[0] = trim_spaces($fields[0]);
    $fields[1] = trim_spaces($fields[1]);
    $metapar{$fields[0]} = $fields[1];
  }
close(META);

my $start_codon_length=$metapar{start_length};
my $start_codon_offset=$metapar{start_offset};
my $intron_short_length = $metapar{intron_short_length};
my $intergenic_length = $metapar{intergenic_length};
my $stop_codon_length= $metapar{stop_codon_length};
my $stop_codon_offset=$metapar{stop_codon_offset};
my $acceptor_length= $metapar{acceptor_length};
my $initial_pattern_length = $metapar{acceptor_offset};
my $acceptor_initial_pattern_length = $metapar{acceptor_initial_pattern_length};
my $donor_initial_pattern_length = $metapar{donor_initial_pattern_length};
my $acceptor_offset=$metapar{acceptor_offset};
my $donor_length=$metapar{donor_length};
my $donor_offset=$metapar{donor_offset};
my $start_initial_pattern_length=$metapar{start_initial_pattern_length};
my $branch_length = $metapar{branch_length};

# Fixing offsets
my $fixed_donor_offset = $donor_offset + $donor_initial_pattern_length;
my $fixed_stop_offset = $stop_codon_offset;
my $donor_signal_length = $donor_length + $donor_initial_pattern_length;
my $acceptor_signal_length = $branch_length + $acceptor_length + $acceptor_initial_pattern_length;
my $start_signal_length = $start_codon_length + 3 + $start_initial_pattern_length;
my $stop_signal_length = $stop_codon_length;
my $exon_length_start = $start_codon_length - $start_codon_offset + $start_initial_pattern_length;
my $exon_length_stop = $stop_codon_offset;
my $exon_length_acceptor =  $acceptor_length - $acceptor_offset - 2 + $acceptor_initial_pattern_length;
my $exon_length_donor = $donor_offset + $donor_initial_pattern_length;
my $exon_delta_initial = $exon_length_start + $exon_length_donor;
my $exon_delta_internal = $exon_length_acceptor + $exon_length_donor;
my $exon_delta_final = $exon_length_acceptor + $exon_length_stop;
my $exon_delta_single = $exon_length_start + $exon_length_stop;
my $intron_delta = $branch_length + $acceptor_offset + 2 + $donor_length- $donor_offset;
my $intron_short_offset_forward = $donor_length - $donor_offset;
my $intron_short_offset_reverse = $branch_length + $acceptor_offset + 2 + $intron_short_length;
my $branch_offset = $branch_length + $acceptor_offset;

# fix the phase of each exon state
my $exon_initial_iphase =  iphase_initial (0);
my $exon_initial_0_ophase =  ophase (0);
my $exon_initial_1_ophase =  ophase (1);
my $exon_initial_2_ophase =  ophase (2);
my $exon_single_ophase = ophase_final(2);
my $exon_final_ophase = ophase_final(2);
my $exon_final_0_iphase = iphase(0);
my $exon_final_1_iphase = iphase(1);
my $exon_final_2_iphase = iphase(2);
my $exon_internal_00_iphase = iphase (0);
my $exon_internal_01_iphase = iphase (0);
my $exon_internal_02_iphase = iphase (0);
my $exon_internal_10_iphase = iphase (1);
my $exon_internal_11_iphase = iphase (1);
my $exon_internal_12_iphase = iphase (1);
my $exon_internal_20_iphase = iphase (2);
my $exon_internal_21_iphase = iphase (2);
my $exon_internal_22_iphase = iphase (2);
my $exon_internal_00_ophase = ophase (0);
my $exon_internal_01_ophase = ophase (1);
my $exon_internal_02_ophase = ophase (2);
my $exon_internal_10_ophase = ophase (0);
my $exon_internal_11_ophase = ophase (1);
my $exon_internal_12_ophase = ophase (2);
my $exon_internal_20_ophase = ophase (0);
my $exon_internal_21_ophase = ophase (1);
my $exon_internal_22_ophase = ophase (2);
my $rexon_initial_ophase =  rophase_initial (0);
my $rexon_initial_0_iphase =  riphase (0);
my $rexon_initial_1_iphase =  riphase (1);
my $rexon_initial_2_iphase =  riphase (2);
my $rexon_single_iphase = riphase_final(2);
my $rexon_single_ophase = rophase_initial(0);
my $rexon_final_0_ophase = rophase(0);
my $rexon_final_1_ophase = rophase(1);
my $rexon_final_2_ophase = rophase(2);
my $rexon_final_iphase = riphase_final(2);
my $rexon_internal_00_iphase = riphase (0);
my $rexon_internal_01_iphase = riphase (0);
my $rexon_internal_02_iphase = riphase (0);
my $rexon_internal_10_iphase = riphase (1);
my $rexon_internal_11_iphase = riphase (1);
my $rexon_internal_12_iphase = riphase (1);
my $rexon_internal_20_iphase = riphase (2);
my $rexon_internal_21_iphase = riphase (2);
my $rexon_internal_22_iphase = riphase (2);
my $rexon_internal_00_ophase = rophase (0);
my $rexon_internal_01_ophase = rophase (1);
my $rexon_internal_02_ophase = rophase (2);
my $rexon_internal_10_ophase = rophase (0);
my $rexon_internal_11_ophase = rophase (1);
my $rexon_internal_12_ophase = rophase (2);
my $rexon_internal_20_ophase = rophase (0);
my $rexon_internal_21_ophase = rophase (1);
my $rexon_internal_22_ophase = rophase (2);


my $ghmm =<<EOF;
model_name = "GeneralizedHiddenMarkovModel"
state_names = ("N","I0","I1","I2",
"EI0","EI1","EI2","ES",
"ET0","ET1","ET2",
"start", "stop", "acc0", "acc1", "acc2", "don0", "don1", "don2",
"E00", "E01", "E02", "E10", "E11", "E12","E20", "E21", "E22",
"rI0","rI1", "rI2",
"rEI0","rEI1","rEI2","rES",
"rET0","rET1","rET2",
"rstart", "rstop", "racc0", "racc1", "racc2", "rdon0", "rdon1", "rdon2",
"rE00", "rE01", "rE02", "rE10", "rE11", "rE12","rE20", "rE21", "rE22", "F"
)
observation_symbols = ("A", "C", "G", "T")
initial_probabilities = ("N": 1.0)
terminal_probabilities = ("F": 1.0)

transitions = (
)

acceptor_model =  "../ghmm/model/acceptor_composed.model"
racceptor_model =  "../ghmm/model/acceptor_composed_rev.model"
donor_model = "../ghmm/model/donor_composed.model"
rdonor_model = "../ghmm/model/donor_composed_rev.model"
start_model = "../ghmm/model/start_composed.model"
stop_model = "../ghmm/model/stop.model"
rdonor_model = "../ghmm/model/donor_rev.model"
rstart_model = "../ghmm/model/start_composed_rev.model"
rstop_model = "../ghmm/model/stop_rev.model"
cds_model = "model/cds.model"
rcds_model =  "model/cds_rev.model"
noncoding_model = "model/intergenic.model"
initial_exon_duration = "../ghmm/model/initial.model"
internal_exon_duration = "../ghmm/model/internal.model"
final_exon_duration = "../ghmm/model/final.model"
single_duration = "../ghmm/model/single.model"
intron_length = "../ghmm/model/intron_length.model"

EI_duration_0 = [ model_name="PhasedRunLengthDistribution"
                  input_phase = 0
                  output_phase = 0
                  number_of_phases = 3
                  delta = $exon_delta_initial
                  model = initial_exon_duration ]
EI_duration_1 = [ model_name="PhasedRunLengthDistribution"
                input_phase = 0
                output_phase = 1
                number_of_phases = 3
                delta = $exon_delta_initial
                model =initial_exon_duration]
EI_duration_2 = [ model_name="PhasedRunLengthDistribution"
                input_phase = 0
                output_phase = 2
                number_of_phases = 3
                delta = $exon_delta_initial
                model = initial_exon_duration]
ES_duration = [ model_name="PhasedRunLengthDistribution"
                 input_phase = 0
                output_phase = 2
                number_of_phases = 3
                delta = $exon_delta_single
                model = single_duration ]
ET0_duration = [ model_name="PhasedRunLengthDistribution"
                input_phase = 0
                output_phase = 2
                number_of_phases = 3
                delta = $exon_delta_final
                model = final_exon_duration ]
ET1_duration = [ model_name="PhasedRunLengthDistribution"
                input_phase = 1
                output_phase = 2
                number_of_phases = 3
                delta =  $exon_delta_final
               model = final_exon_duration]
ET2_duration = [ model_name="PhasedRunLengthDistribution"
                input_phase = 2
                output_phase = 2
                number_of_phases = 3
                delta = $exon_delta_final
                model =final_exon_duration ]
E0_duration = [ model_name="PhasedRunLengthDistribution"
                input_phase = 0
                output_phase = 0
                number_of_phases = 3
                delta = $exon_delta_internal
                model =internal_exon_duration ]
E1_duration = [ model_name="PhasedRunLengthDistribution"
                input_phase = 0
                output_phase = 1
                number_of_phases = 3
                delta = $exon_delta_internal
                model = internal_exon_duration ]
E2_duration = [ model_name="PhasedRunLengthDistribution"
                input_phase = 0
                output_phase = 2
                number_of_phases = 3
                delta = $exon_delta_internal
                model = internal_exon_duration ]
EI0 = [ observation = cds_model
        duration = EI_duration_0
        extend_emission = ($exon_length_start, $exon_length_donor)
        input_phase = $exon_initial_iphase
        output_phase = $exon_initial_0_ophase]
EI1 = [ observation = cds_model
        duration = EI_duration_1
        extend_emission = ($exon_length_start, $exon_length_donor)
        input_phase = $exon_initial_iphase
        output_phase = $exon_initial_1_ophase]
EI2 = [ observation = cds_model
        duration = EI_duration_2
        extend_emission = ($exon_length_start, $exon_length_donor)
        input_phase = $exon_initial_iphase
        output_phase = $exon_initial_2_ophase]
ES = [ observation = cds_model
       duration = ES_duration
       extend_emission = ($exon_length_start, $exon_length_stop)
       input_phase = $exon_initial_iphase
       output_phase = $exon_single_ophase]
ET0 = [ observation = cds_model
        duration = ET0_duration
        extend_emission = ($exon_length_acceptor, $exon_length_stop)
        input_phase = $exon_final_0_iphase
        output_phase = $exon_final_ophase]
ET1 = [ observation = cds_model
        duration = ET1_duration
        extend_emission = ($exon_length_acceptor, $exon_length_stop)
        input_phase = $exon_final_1_iphase
        output_phase = $exon_final_ophase]
ET2 = [ observation  = cds_model
        duration = ET2_duration
        extend_emission = ($exon_length_acceptor, $exon_length_stop)
        input_phase  = $exon_final_2_iphase
        output_phase =  $exon_final_ophase]
E00 = [ observation  = cds_model
        duration = E0_duration
        extend_emission = ($exon_length_acceptor, $exon_length_donor)
        input_phase  = $exon_internal_00_iphase
        output_phase = $exon_internal_00_ophase ]
E01= [ observation  = cds_model
       duration = E1_duration
        extend_emission = ($exon_length_acceptor, $exon_length_donor)
       input_phase  = $exon_internal_01_iphase
       output_phase =  $exon_internal_01_ophase ]
E02= [ observation  = cds_model
       duration = E2_duration
        extend_emission = ($exon_length_acceptor, $exon_length_donor)
       input_phase  = $exon_internal_02_iphase
       output_phase =  $exon_internal_02_ophase ]
E10= [ observation  = cds_model
        duration = E2_duration
        extend_emission = ($exon_length_acceptor, $exon_length_donor)
       input_phase  = $exon_internal_10_iphase
       output_phase =  $exon_internal_10_ophase ]
E11= [ observation  = cds_model
       duration = E0_duration
        extend_emission = ($exon_length_acceptor, $exon_length_donor)
       input_phase  = $exon_internal_11_iphase
       output_phase =  $exon_internal_11_ophase ]
E12= [ observation  = cds_model
        duration = E1_duration
        extend_emission = ($exon_length_acceptor, $exon_length_donor)
        input_phase  = $exon_internal_12_iphase
        output_phase =  $exon_internal_12_ophase ]
E20= [ observation  = cds_model
        duration = E1_duration
        extend_emission = ($exon_length_acceptor, $exon_length_donor)
        input_phase  = $exon_internal_20_iphase
        output_phase =  $exon_internal_20_ophase ]
E21= [ observation  = cds_model
        duration = E2_duration
       extend_emission = ($exon_length_acceptor, $exon_length_donor)
        input_phase  = $exon_internal_21_iphase
        output_phase =  $exon_internal_21_ophase ]
E22= [ observation  = cds_model
       duration = E0_duration
       extend_emission = ($exon_length_acceptor, $exon_length_donor)
       input_phase  = $exon_internal_22_iphase
       output_phase =  $exon_internal_22_ophase ]
start = [ observation = start_model
          sequence_length = $start_signal_length ]
stop = [ observation = stop_model
         sequence_length = $stop_signal_length]
acc0 = [ observation = acceptor_model
          sequence_length = $acceptor_signal_length ]
acc1 = [ observation = acceptor_model
          sequence_length = $acceptor_signal_length ]
acc2 = [ observation = acceptor_model
          sequence_length = $acceptor_signal_length ]
don0 = [ observation = donor_model
          sequence_length = $donor_signal_length ]
don1 = [ observation = donor_model
          sequence_length = $donor_signal_length ]
don2 = [ observation = donor_model
          sequence_length = $donor_signal_length ]
N = [observation = noncoding_model]
F = [observation = noncoding_model]
I0 = [ observation = noncoding_model ]
I1 = [ observation = noncoding_model ]
I2 = [ observation = noncoding_model ]
rI0 = [ observation = noncoding_model ]
rI1 = [ observation = noncoding_model ]
rI2 = [ observation = noncoding_model ]
rEI0  = [ observation = rcds_model
          duration = EI_duration_0
          extend_emission = ($exon_length_donor, $exon_length_start)
          input_phase = $rexon_initial_0_iphase
          output_phase = $rexon_initial_ophase ]
rEI1 = [ observation = rcds_model
         duration = EI_duration_1
         extend_emission = ($exon_length_donor, $exon_length_start)
         input_phase = $rexon_initial_1_iphase
         output_phase = $rexon_initial_ophase ]
rEI2 = [ observation = rcds_model
         duration = EI_duration_2
         extend_emission = ($exon_length_donor, $exon_length_start)
         input_phase = $rexon_initial_2_iphase
         output_phase = $rexon_initial_ophase ]
rES = [ observation =rcds_model
        duration = ES_duration
        extend_emission = ($exon_length_stop, $exon_length_start)
        input_phase = $rexon_single_iphase
        output_phase = $rexon_single_ophase ]
rET0 = [ observation = rcds_model
         duration = ET0_duration
         extend_emission = ($exon_length_stop, $exon_length_acceptor)
         input_phase = $rexon_final_iphase
         output_phase = $rexon_final_0_ophase ]
rET1 = [ observation = rcds_model
         duration = ET1_duration
         extend_emission = ($exon_length_stop, $exon_length_acceptor)
         input_phase = $rexon_final_iphase
         output_phase = $rexon_final_1_ophase ]
rET2 = [ observation  = rcds_model
         duration = ET2_duration
         extend_emission = ($exon_length_stop, $exon_length_acceptor)
         input_phase  = $rexon_final_iphase
         output_phase = $rexon_final_2_ophase ]
rE00 = [ observation  = rcds_model
         duration = E0_duration
         extend_emission = ($exon_length_donor, $exon_length_acceptor)
         input_phase  = $rexon_internal_00_iphase
         output_phase = $rexon_internal_00_ophase ]
rE01= [ observation  = rcds_model
        duration = E2_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_01_iphase
        output_phase = $rexon_internal_01_ophase ]
rE02= [ observation  = rcds_model
        duration = E1_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_02_iphase
        output_phase = $rexon_internal_02_ophase ]
rE10= [ observation  = rcds_model
        duration = E1_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_10_iphase
        output_phase = $rexon_internal_10_ophase ]
rE11= [ observation  = rcds_model
        duration = E0_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_11_iphase
        output_phase = $rexon_internal_11_ophase ]
rE12= [ observation  = rcds_model
        duration = E2_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_12_iphase
        output_phase = $rexon_internal_12_ophase ]
rE20= [ observation  = rcds_model
        duration = E2_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_20_iphase
        output_phase = $rexon_internal_20_ophase ]
rE21= [ observation  = rcds_model
        duration = E1_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_21_iphase
        output_phase = $rexon_internal_21_ophase ]
rE22= [ observation  = rcds_model
        duration = E0_duration
        extend_emission = ($exon_length_donor, $exon_length_acceptor)
        input_phase  = $rexon_internal_22_iphase
        output_phase = $rexon_internal_22_ophase ]
rstart = [ observation = rstart_model
          sequence_length = $start_signal_length ]
rstop = [ observation = rstop_model
          sequence_length = $stop_signal_length ]
racc0 = [ observation = racceptor_model
          sequence_length = $acceptor_signal_length  ]
racc1 = [ observation = racceptor_model
          sequence_length = $acceptor_signal_length  ]
racc2 = [ observation = racceptor_model
          sequence_length = $acceptor_signal_length  ]
rdon0 = [ observation = rdonor_model
          sequence_length = $donor_signal_length ]
rdon1 = [ observation = rdonor_model
          sequence_length = $donor_signal_length ]
rdon2 = [ observation = rdonor_model
          sequence_length = $donor_signal_length ]
EOF
open (GHMM, ">ghmm/model/ghmm_geometric_nostop.model") or die "model/ghmm_geometric_nostop.model: $!";
print GHMM $ghmm;
close(GHMM);

###################

my $final_state_prob = $metapar{final_state_probability};


my %estimated;
open (MODEL, "< ghmm/model/states.model" ) ;
foreach my $line (<MODEL>) {
    if($line =~ m/\|/) {
        $line =~ s/\s//g;
        $line =~ s/;//g;
        $line =~ s/#leaf//g;
        if(!($line  =~ m/0;/))
        {
            my ($trans, $probs) = split(":", $line);
            $estimated{$trans} = $probs;
        }
    }
}
close(MODEL);

my $pforward = 0.5;

my $intron_lengths =`myop-sequence_length  < ghmm/dataset/intron.fasta`;
$intron_lengths =~ s/(.+):\t//g;
my @intron_lengths = split(/\s/, $intron_lengths);
my $count = 0;
my $length_mean= 0;
my @ilengths;

my $num_short_introns = 0;
my $num_long_introns = 0;
my %counter;
my $max_intron_length = -10;
my $sum_intron_length = 0;
foreach my $l (@intron_lengths)
 {
    $counter{$l} ++;
    if($l > 0) {
      $sum_intron_length += $l;
    }
    if($max_intron_length < $l)
    {
        $max_intron_length = $l;
    }
}

my $total_introns = scalar(@intron_lengths);
for(my $i = 0; $i <= $intron_short_length; $i++)
{
    if(defined $counter{$i}) {
        $num_short_introns += $counter{$i};
    }
}

for(my $i = $intron_short_length+1; $i <= $max_intron_length; $i++)
{
    if(defined $counter{$i}) {
        $num_long_introns += $counter{$i} * $i;
    }
}


my $seq = "seq: $intron_short_length\n";
my ($prob_length_d, @lixo) = split("\n", `echo "$seq" | tops-evaluate -m ghmm/model/intron_length.model | awk -F" " '{print\$2}'`);
my $probd = exp($prob_length_d);

print STDERR "P(m=d) = ".$probd."\n";
my $mlong = $num_long_introns/ ($total_introns - $num_short_introns) - $intron_short_length;
if (($mlong < 0) || ($total_introns - $num_short_introns) == 0) {
    $mlong = 10;
}

my $mean_intron_length = ($sum_intron_length/scalar(@intron_lengths)) ;
my @lengths = sort {$a <=> $b} (@intron_lengths);

my $pshort = (1.0)/ (1.0 + $probd * $mlong);
my $p = 1.0 - 1.0/$mlong;

my $pd = $probd * $pshort;
my $pd1 = (1.0 - $pshort)/$mlong;

print STDERR "freq_short = ".((($num_short_introns)/scalar(@intron_lengths))*100.0)."\n";
print STDERR "pshort = $pshort\n";
print STDERR "q = ".(1-$p)."\n";
print STDERR "mlong = $mlong\n";
print STDERR "P(L=d) = $pd\n";
print STDERR "P(L=d+1) = $pd1\n";

my %transitions;


$transitions{make_transition_str("don0", "I0")} = 1.0;
$transitions{make_transition_str("I0", "I0")} = $p;
$transitions{make_transition_str("I0", "acc0")} = 1.0 - $p;

$transitions{make_transition_str("don1", "I1")} = 1.0;
$transitions{make_transition_str("I1", "I1")} = $p;
$transitions{make_transition_str("I1", "acc1")} = 1.0 - $p;

$transitions{make_transition_str("don2", "I2")} = 1.0;
$transitions{make_transition_str("I2", "I2")} = $p;
$transitions{make_transition_str("I2", "acc2")} = 1.0 - $p;

$transitions{make_transition_str("racc0", "rI0")} = 1.0;
$transitions{make_transition_str("rI0", "rI0")} = $p;
$transitions{make_transition_str("rI0", "rdon0")} = 1.0 - $p;
$transitions{make_transition_str("racc1", "rI1")} = 1.0;
$transitions{make_transition_str("rI1", "rI1")} = $p;
$transitions{make_transition_str("rI1", "rdon1")} = 1.0 - $p;
$transitions{make_transition_str("racc2", "rI2")} = 1.0;
$transitions{make_transition_str("rI2", "rI2")} = $p;
$transitions{make_transition_str("rI2", "rdon2")} = 1.0 - $p;
$transitions{make_transition_str( "F", "F" )} =  1.0;
$transitions{make_transition_str( "N", "start" )} = $pforward * ((1.0/$intergenic_length) - $final_state_prob);
$transitions{make_transition_str( "N", "rstop" )} = (1.0 - $pforward) * ((1.0/$intergenic_length) - $final_state_prob);
$transitions{make_transition_str("N", "N")} = 1.0 - (1.0/$intergenic_length);
$transitions{make_transition_str( "N", "F" )} = $final_state_prob;
$transitions{make_transition_str( "stop", "N" )} = 1;
$transitions{make_transition_str( "rstart", "N" ) }= 1;
$transitions{make_transition_str("start","ES")} = 0.1;
$transitions{make_transition_str("ES","stop")} = 1.0;
$transitions{make_transition_str("rstop","rES")} = 0.1;
$transitions{make_transition_str("rES","rstart")} = 1.0;
$transitions{make_transition_str("start", "EI0")} = 0.3;
$transitions{make_transition_str("start", "EI1")} = 0.3;
$transitions{make_transition_str("start", "EI2")} = 0.3;
$transitions{make_transition_str("EI0", "don1")} = 1.0;
$transitions{make_transition_str("EI1", "don2")} = 1.0;
$transitions{make_transition_str("EI2", "don0")} = 1.0;
$transitions{make_transition_str("acc0", "ET0")} = 0.1;
$transitions{make_transition_str("acc0", "E00")} = 0.3;
$transitions{make_transition_str("acc0", "E01")} = 0.3;
$transitions{make_transition_str("acc0", "E02")} = 0.3;
$transitions{make_transition_str("acc1", "ET1")} = 0.1;
$transitions{make_transition_str("acc1", "E10")} = 0.3;
$transitions{make_transition_str("acc1", "E11")} = 0.3;
$transitions{make_transition_str("acc1", "E12")} = 0.3;
$transitions{make_transition_str("acc2", "ET2")} = 0.1;
$transitions{make_transition_str("acc2", "E20")} = 0.3;
$transitions{make_transition_str("acc2", "E21")} = 0.3;
$transitions{make_transition_str("acc2", "E22")} = 0.3;
$transitions{make_transition_str("ET0", "stop")} = 1.0;
$transitions{make_transition_str("ET1", "stop")} = 1.0;
$transitions{make_transition_str("ET2", "stop")} = 1.0;
$transitions{make_transition_str("E00", "don1")} = 1.0;
$transitions{make_transition_str("E01", "don2")} = 1.0;
$transitions{make_transition_str("E02", "don0")} = 1.0;
$transitions{make_transition_str("E10", "don1")} = 1.0;
$transitions{make_transition_str("E11", "don2")} = 1.0;
$transitions{make_transition_str("E12", "don0")} = 1.0;
$transitions{make_transition_str("E20", "don1")} = 1.0;
$transitions{make_transition_str("E21", "don2")} = 1.0;
$transitions{make_transition_str("E22", "don0")} = 1.0;
$transitions{make_transition_str("rstop","rET0")} = 0.3;
$transitions{make_transition_str("rstop","rET1")} = 0.3;
$transitions{make_transition_str("rstop","rET2")} = 0.3;
$transitions{make_transition_str("rstop","rES")} = 0.1;
$transitions{make_transition_str("rET2","racc2")} = 1.0;
$transitions{make_transition_str("rET0","racc0")} = 1.0;
$transitions{make_transition_str("rET1","racc1")} = 1.0;
$transitions{make_transition_str("rdon1","rEI0")} = 0.1;
$transitions{make_transition_str("rdon1","rE00")} = 0.3;
$transitions{make_transition_str("rdon1","rE01")} = 0.3;
$transitions{make_transition_str("rdon1","rE02")} = 0.3;
$transitions{make_transition_str("rdon0","rEI2")} = 0.1;
$transitions{make_transition_str("rdon0","rE20")} = 0.3;
$transitions{make_transition_str("rdon0","rE21")} = 0.3;
$transitions{make_transition_str("rdon0","rE22")} = 0.3;
$transitions{make_transition_str("rdon2","rEI1")} = 0.1;
$transitions{make_transition_str("rdon2","rE10")} = 0.3;
$transitions{make_transition_str("rdon2","rE11")} = 0.3;
$transitions{make_transition_str("rdon2","rE12")} = 0.3;
$transitions{make_transition_str("rE00","racc0")} = 1.0;
$transitions{make_transition_str("rE01","racc1")} = 1.0;
$transitions{make_transition_str("rE02","racc2")} = 1.0;
$transitions{make_transition_str("rE10","racc0")} = 1.0;
$transitions{make_transition_str("rE11","racc1")} = 1.0;
$transitions{make_transition_str("rE12","racc2")} = 1.0;
$transitions{make_transition_str("rE20","racc0")} = 1.0;
$transitions{make_transition_str("rE21","racc1")} = 1.0;
$transitions{make_transition_str("rE22","racc2")} = 1.0;
$transitions{make_transition_str("rEI0","rstart")} = 1.0;
$transitions{make_transition_str("rEI1","rstart")} = 1.0;
$transitions{make_transition_str("rEI2","rstart")} = 1.0;



open (GHMM2, ">ghmm/model/ghmm_final.model") or die "$!";
open (GHMM, "<ghmm/model/ghmm_geometric_nostop.model") or die "$!";
my $reading_trans = 0;
foreach my $line (<GHMM>) {
  if ($line =~ m/transitions =/) {
    print GHMM2 "transitions = (\n";
    foreach my $trans (sort
                       {
                         my ($a1, $a2) = split(/\|/, $a);
                         my ($b1, $b2) = split(/\|/, $b);
                         return $a2.$a1 cmp $b2.$b1;

                       } (keys %transitions)  ) {

      print GHMM2 $trans.": ".$transitions{$trans}.";\n";
    }
    print GHMM2 ")\n";
    $reading_trans = 1;
  }
  if(!$reading_trans) {
    print GHMM2 $line;
  }
  if($reading_trans && $line =~ /\)/) {
    $reading_trans = 0;
  }

}
close(GHMM2);
close(GHMM);

system("rm ghmm/model/ghmm_geometric_nostop.model");
system("mv ghmm/model/ghmm_final.model ghmm/model/ghmm_geometric_nostop.model");


sub make_transition_str {
    my $from = shift;
    my $to = shift;
    return "\"".$to."\"|\"".$from."\"";
}



sub state_sequence_forward {

    my $cds_lengths = shift;
    my $reversed = shift;
    my @cds = @$cds_lengths;
    my $state_seq = " N";
    my $i = 0;
    my $k = 1;
    if($reversed) {
        $i = $#cds;
        $k = -1;
    }
    my $iphase = 0;
    my $ophase = 0;
    if(scalar (@cds) > 1 ) {
        $state_seq .= " start EI";
        my $length = $cds[$i];
        $ophase = ($length + $iphase - 1) % 3;
        $iphase = ($ophase + 1) % 3;
        $state_seq .= $ophase;
        $state_seq .= " don$iphase I$iphase acc$iphase";
        $i += $k;
        while(($i < $#cds) && ($i > 0))
        {
            $state_seq .= " E$iphase";
            my $length = $cds[$i];
            $ophase = ($length + $iphase -1) % 3;
            $iphase = ($ophase + 1) % 3;
            $state_seq .= $ophase;
            $state_seq .= " don$iphase I$iphase acc$iphase";
            $i += $k;
        }
        $state_seq .= " ET$iphase";
        $length = $cds[$i];
        $ophase = ($length + $iphase - 1) % 3;
        $iphase = ($ophase + 1) % 3;
        $state_seq .= " stop";
    } elsif ($#cds == 0) {
        $state_seq .= " start ES stop N";
    }

    return $state_seq;

}

sub state_sequence_reverse {
    my $cds_lengths = shift;
    my $reversed = shift;
    my @cds = @$cds_lengths;
    my $iphase = 0;
    my $ophase = 0;
    my $state_seq  = " N";
    my $i = 0;
    my $k = 1;
    if($reversed) {
        $i = $#cds;
        $k = -1;
    }
    if(scalar (@cds) > 1 ) {
        $state_seq .= " rstop rET";
        my $length = $cds[$i];
        $ophase = (3-($length + $iphase-1 )%3)%3;
        $iphase = ($ophase - 1) % 3;
        $state_seq .= $ophase;
        $state_seq .= " racc$iphase rI$iphase rdon$iphase";
        $i += $k;
        while(($i < $#cds) && ($i > 0) )
        {
            $state_seq .= " rE$iphase";
            my $length = $cds[$i];
            $ophase = (3 - ($length + $iphase -1)%3) % 3;
            $iphase = ($ophase - 1) % 3;
            $state_seq .= $ophase;
            $state_seq .= " racc$iphase rI$iphase rdon$iphase";
            $i += $k;
        }
        $state_seq .= " rEI$iphase";
        $length = $cds[$i];
        $ophase = (3 - ($length + $iphase -1) % 3) % 3;
        $iphase = ($ophase - 1) % 3;
        $state_seq .= " rstart";
    } elsif ($#cds == 0) {
        $state_seq .= " rstop rES rstart N";
    }
    return $state_seq;
}




####################

sub iphase_initial {
  return $exon_length_start % 3;
}

sub iphase {
  my $iphase = shift;
  return ($exon_length_acceptor + $iphase) % 3;
}

sub ophase {
  my $ophase = shift;
  return ($ophase - $exon_length_donor) % 3;
}

sub ophase_final {
  my $ophase = shift;
  return ($ophase - $exon_length_stop ) % 3;
}

sub rophase_initial {
  my $ophase = shift;
  my $phase = 1;
  if($ophase == 2) {
    $phase = 0;
  } elsif($ophase == 0){
    $phase = 2;
  }
  return ($phase - $exon_length_start) % 3;
}


sub riphase {
  my $iphase = shift;
  my $phase = 1;
  if($iphase == 2) {
    $phase = 0;
  } elsif($iphase == 0){
    $phase = 2;
  }
  return ($exon_length_donor + $phase) % 3;
}

sub rophase {
  my $ophase = shift;
  my $phase = 1;
  if($ophase == 2) {
    $phase = 0;
  } elsif($ophase == 0){
    $phase = 2;
  }
  return ($phase - $exon_length_acceptor ) % 3;
}

sub riphase_final {
  my $ophase = shift;
  my $phase = 1;
  if($ophase == 2) {
    $phase = 0;
  } elsif($ophase == 0){
    $phase = 2;
  }
  return ($exon_length_stop + $phase) % 3;
}


sub trim_spaces {
  my $v = shift;
  $v =~ s/^\s+//;     $v =~ s/\s+$//;
  return $v;
}
