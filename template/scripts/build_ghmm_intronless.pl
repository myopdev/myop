#!/usr/bin/perl

open (OUT, ">ghmm/model/ghmm_intronless.model") or die "model/ghmm_intronless.model: $!\n";

print OUT <<EOF;

model_name = "GeneralizedHiddenMarkovModel"
state_names = ("begin", "end", "CDS",   "start",  "stop", "CDS0", "CDS1", "CDS2")
observation_symbols = ("A", "C", "G", "T")
initial_probabilities = ("begin": 0.999; "CDS2": 0.0003; "CDS1": 0.0003; "CDS0": 0.0004;  )

transitions = ( 
"begin" | "begin": 0.9999;
"start" | "begin": 0.0001;
"CDS" | "start": 1;
"end" | "stop": 1;
"stop" | "CDS": 1;
"end" | "end": 1;
"stop" | "CDS0": 1;
"stop" | "CDS1": 1;
"stop" | "CDS2": 1;
)
start_model = "../ghmm/model/start_composed_motifless.model"
stop_model = "../ghmm/model/stop.model"
cds_model = "model/cds.model"
noncoding_model = "model/intergenic.model"
cds_length = "../ghmm/model/cds_length.model"

cds_duration = [ model_name="PhasedRunLengthDistribution"
                 input_phase = 0
                 output_phase = 2
                 number_of_phases = 3
                 delta = 4
                 model = cds_length ]
cds_duration_0 = [ model_name="PhasedRunLengthDistribution"
                 input_phase = 0
                 output_phase = 2
                 number_of_phases = 3
                 delta = 0
                 model = cds_length ]

cds_duration_1 = [ model_name="PhasedRunLengthDistribution"
                 input_phase = 1
                 output_phase = 2
                 number_of_phases = 3
                 delta = 0
                 model = cds_length ]
cds_duration_2 = [ model_name="PhasedRunLengthDistribution"
                 input_phase = 2
                 output_phase = 2
                 number_of_phases = 3
                 delta = 0
                 model = cds_length ]


start = [ observation = start_model
          sequence_length = 7 ]
stop = [ observation = stop_model
         sequence_length = 3]
begin = [ observation = noncoding_model ]
end = [ observation = noncoding_model ]
CDS0 = [ observation = cds_model 
         duration = cds_duration_0
         input_phase =  0 
         output_phase = 2] 
CDS1 = [ observation = cds_model 
         duration = cds_duration_1
         input_phase =  1
         output_phase = 2] 
CDS2 = [ observation = cds_model 
         duration = cds_duration_2
         input_phase =  2
         output_phase = 2] 
CDS = [ observation = cds_model 
        duration = cds_duration
        extend_emission = (4, 0)
        input_phase =  1
        output_phase = 2 ] 
EOF
close(OUT);
