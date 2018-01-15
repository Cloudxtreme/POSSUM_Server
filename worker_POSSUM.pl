#!/usr/bin/perl
#use strict;
#use warnings;
use Gearman::Worker;
use Storable qw(thaw);
use Storable qw(freeze);
use DBI;
use IO::All;
use Bio::SeqIO;
use Bio::Seq;
use Capture::Tiny ':all';

print "start new worker \n";
my $worker = Gearman::Worker->new;
print "finish new worker \n";

print "start job_servers \n";
$worker->job_servers('127.0.0.1',4730);
print "finish job_servers \n";

print "start register_function \n";
$worker->register_function( generatePSSM => \&generatePSSM );
print "finish register_function \n";

print "start work \n";
$worker->work while 1;
print "finish work \n";


sub generatePSSM
{
    print "enter the function generatePSSM \n";
    my @row=@{ thaw($_[0]->arg) };

    print "start print the args \n";
    print "\n";
    print "$row[0] \n";
    print "$row[1] \n";
    print "$row[2] \n";
    print "$row[3] \n";
    print "$row[4] \n";
    print "$row[5] \n";
    print "$row[6] \n";
    print "$row[7] \n";
    print "$row[8] \n";
    print "$row[9] \n";
    print "\n";
    print "finish print the args \n";

    my $class = $row[0];
    my $jobName = $row[1];
    my $seqs = $row[2];
    my $seqNum = $row[3];
    my $selectedFeatures = $row[4];
    my $email = $row[5];
    my $databaseType = $row[6];
    my $arguments = $row[7];
    my $iterations = $row[8];
    my $E_value = $row[9];

    my $blastpgp = "/possum/softwares/blast-2.2.26/bin/blastpgp";
    my $blastpgp_db = "";
    my $pssm_cache_table = "";
    my $pssm_cache_dir = "";
    my $tomcat_storage_dir = "/possum/tomcatWebapps/ROOT/static/result_files";
    if($databaseType==1)
    {
        $pssm_cache_table = "pssm_uniref50";
        $pssm_cache_dir = "/possum/POSSUM_cache/pssm_data_uniref50";
        $blastpgp_db = "/possum/uniref/uniref50/uniref50";
    }
    if($databaseType==2)
    {
        $pssm_cache_table = "pssm_uniref90";
        $pssm_cache_dir = "/possum/POSSUM_cache/pssm_data_uniref90";
        $blastpgp_db = "/possum/uniref/uniref90/uniref90.fasta";
    }
    if($databaseType==3)
    {
        $pssm_cache_table = "pssm_uniref100";
        $pssm_cache_dir = "/possum/POSSUM_cache/pssm_data_uniref100";
        $blastpgp_db = "/possum/uniref/uniref100/uniref100.fasta";
    }


    #################################updateStatusAndStartTimeInJob
    print "worker_POSSUM -- begin to update status in table job:\n";
    my $STATUS_PROCESSING = 1;
    my $START_TIME = 1;
    &updateStatus($jobName,$STATUS_PROCESSING,$START_TIME);
    print "worker_POSSUM -- Success to update status in table job!\n";
    #################################updateStatusAndStartTimeInJob

    local $retStr="";
    local $errorMsg="";
    my $result_folder = "$tomcat_storage_dir/$jobName";
    if (-e "$result_folder")
    {
        print "$result_folder has existed!\n";
        print "sudo rm -fr $result_folder\n";
        `sudo rm -fr $result_folder`;
    }
    print "sudo mkdir $result_folder\n";
    `sudo mkdir $result_folder`;

    my @featureVec = split(//,$selectedFeatures);
    my @argsVec= split(/;/,$arguments);
    print "Begin to store sequences into sequence/$jobName.fasta\n";
    $seqs > io("sequence/$jobName.fasta");
    print "Success to store sequences into sequence/$jobName.fasta\n";

    my $folder = "seperatedSequence/$jobName";
    if (-e "$folder")
    {
        print "$folder has existed!\n";
        print "sudo rm -fr $folder\n";
        `sudo rm -fr $folder`;
    }
    print "sudo mkdir $folder\n";
    `sudo mkdir $folder`;

    print "Begin to divide sequence/$jobName.fasta into seperatedSequence/$jobName:\n";
    $errorMsg = tee_merged {
        `sudo featureCalcCode/divideFastaIntoSingleEntry.pl sequence/$jobName.fasta seperatedSequence/$jobName $jobName`;
    };

    if( ( $? >> 8 ) != 0 ) {
        $retStr = $retStr."error;";
        $retStr = $retStr.$jobName;
        $retStr = $retStr.";";
        $retStr = $retStr.$selectedFeatures;
        $retStr = $retStr.";";
        $retStr = $retStr.$errorMsg;
	print "$retStr";
        print "\n";
        return $retStr;
    }
    print "Success to divide sequence/$jobName.fasta into seperatedSequence/$jobName\n";

    my $newSequences = "";
    my $sequences_missed_pssm = "";
    my $missed_pssm_count = 0;

    print "Begin the workflow\n";
    my $beginWorkflowTime = time();
    ###########################################################
    ################   Begin to generate PSSM profiles  ####################
    ###########################################################
    print "Begin to generate pssm file\n";
    my $beginGenerateTime = time();

    my $folder_in_pssm = "pssm/$jobName";
    if (-e "$folder_in_pssm")
    {
        print "$folder_in_pssm has existed!\n";
        print "sudo rm -fr $folder_in_pssm\n";
        `sudo rm -fr $folder_in_pssm`;
    }
    print "sudo mkdir $folder_in_pssm\n";
    `sudo mkdir $folder_in_pssm`;


    foreach (1..$seqNum)
    {
        ##############################################
        my $seq = &getSeqByPath("seperatedSequence/$jobName/$jobName\_$_.fasta");
        my $pssmName = &getCacheFileName($pssm_cache_table,$seq);
        if(($iterations==3)&&($E_value==0.001)&&$pssmName)
        {
            print "The pssm file for $jobName\_$_.fasta has exsited!\n";
            print "sudo cp $pssm_cache_dir/$pssmName pssm/$jobName/$jobName\_$_.pssm\n";

            $errorMsg = tee_merged {
                `sudo cp $pssm_cache_dir/$pssmName pssm/$jobName/$jobName\_$_.pssm`;
            };

            if( ( $? >> 8 ) != 0 ) {
                $retStr = $retStr."error;";
                $retStr = $retStr.$jobName;
                $retStr = $retStr.";";
                $retStr = $retStr.$selectedFeatures;
                $retStr = $retStr.";";
                $retStr = $retStr.$errorMsg;
                return $retStr;
            }
        }
        else
        {
            print "Begin to generate pssm file for $jobName\_$_.fasta:\n";
            #########################################################
            ################### Need modify   #######################
            #########################################################
            print "sudo $blastpgp -i seperatedSequence/$jobName/$jobName\_$_.fasta -d $blastpgp_db -h $E_value -e 0.0001 -j $iterations -Q pssm/$jobName/$jobName\_$_.pssm\n";

            $errorMsg = tee_merged {
                `sudo $blastpgp -i seperatedSequence/$jobName/$jobName\_$_.fasta -d $blastpgp_db -h $E_value -e 0.0001 -j $iterations -Q pssm/$jobName/$jobName\_$_.pssm`;
            };

            if( ( $? >> 8 ) != 0 ) {
                $retStr = $retStr."error;";
                $retStr = $retStr.$jobName;
                $retStr = $retStr.";";
                $retStr = $retStr.$selectedFeatures;
                $retStr = $retStr.";";
                $retStr = $retStr.$errorMsg;
                return $retStr;
            }
            print "Finish to generate pssm file for $jobName\_$_.fasta:\n";
            print "sudo cp pssm/$jobName/$jobName\_$_.pssm $pssm_cache_dir/\n";

            $errorMsg = tee_merged {
                `sudo cp pssm/$jobName/$jobName\_$_.pssm $pssm_cache_dir/`;
            };

            if( ( $? >> 8 ) != 0 ) {
                $retStr = $retStr."error;";
                $retStr = $retStr.$jobName;
                $retStr = $retStr.";";
                $retStr = $retStr.$selectedFeatures;
                $retStr = $retStr.";";
                $retStr = $retStr.$errorMsg;
                return $retStr;
            }

            if(($iterations==3)&&($E_value==0.001)){
            	print "worker_POSSUM -- Begin to record the pssm file name into database:\n";
            	&storeCacheFileName($pssm_cache_table,$seq,"$jobName\_$_.pssm");
            	print "worker_POSSUM -- Success to record the pssm file name into database!\n";
            }
        }

        my $output_file = "pssm/$jobName/$jobName\_$_.pssm";
        my $output_file_exist = -e $output_file;
        if(!$output_file_exist)
        {
            my $errorSequence = &getTitleAndSeqByPath("seperatedSequence/$jobName/$jobName\_$_.fasta");
            $sequences_missed_pssm = $sequences_missed_pssm.$errorSequence;
            $missed_pssm_count = $missed_pssm_count+1;
            if(index($retStr,"warning") == -1)
            {
                $retStr = $retStr."warning;";
            }

        }
        ##############################################
    }

    print "Success to generate pssm file\n";
    my $finishGenerateTime = time();

    my $generateInterval = &getMinutes($beginGenerateTime, $finishGenerateTime);
    print "Generating PSSM consuming time is $generateInterval minutes\n";

    ###########################################################
    ################   Finish to generate PSSM profiles  ####################
    ###########################################################

    ##################################################
    ####      Begin to calculate features         ####
    ##################################################

    print "Begin to generate original features\n";
    my $beginCalculatingTime = time();
    my $folder_in_features = "originalFeatures/$jobName";
    if (-e "$folder_in_features")
    {
        print "$folder_in_features has existed!\n";
        print "sudo rm -fr $folder_in_features\n";
        `sudo rm -fr $folder_in_features`;
    }
    print "sudo mkdir $folder_in_features\n";
    `sudo mkdir $folder_in_features`;
    if($featureVec[0]==1)
    {
	    my @args = split(/,/,$argsVec[0]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
		my $now1 = time();
		print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_aac_pssm_no_header.csv -t aac_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";
        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_aac_pssm_no_header.csv -t aac_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }


		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_aac_pssm_no_header.csv -p aac_pssm -n 20 -o originalFeatures/$jobName/$jobName\_aac_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_aac_pssm_no_header.csv -p aac_pssm -n 20 -o originalFeatures/$jobName/$jobName\_aac_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }

  		print "sudo rm originalFeatures/$jobName/$jobName\_aac_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_aac_pssm_no_header.csv`;
		my $now2 = time();
		my $minutesInterval = &getMinutes($now1, $now2);
		print "AAC-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[1]==1)
    {
	    my @args = split(/,/,$argsVec[1]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
       	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_d_fpssm_no_header.csv -t d_fpssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_d_fpssm_no_header.csv -t d_fpssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }

       	print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_d_fpssm_no_header.csv -p d_fpssm -n 20 -o originalFeatures/$jobName/$jobName\_d_fpssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_d_fpssm_no_header.csv -p d_fpssm -n 20 -o originalFeatures/$jobName/$jobName\_d_fpssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }

  		print "sudo rm originalFeatures/$jobName/$jobName\_d_fpssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_d_fpssm_no_header.csv`;
      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "D-FPSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[2]==1)
    {
	    my @args = split(/,/,$argsVec[2]);
    	my $arg0 = $args[0];
    	my $arg1 = $args[1];
      	my $now1 = time();
      	my $dimension2 = $arg1*20;
       	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_smoothed_pssm_no_header.csv -t smoothed_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_smoothed_pssm_no_header.csv -t smoothed_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }

       	print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_smoothed_pssm_no_header.csv -p smoothed_pssm -n $dimension2 -o originalFeatures/$jobName/$jobName\_smoothed_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_smoothed_pssm_no_header.csv -p smoothed_pssm -n $dimension2 -o originalFeatures/$jobName/$jobName\_smoothed_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_smoothed_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_smoothed_pssm_no_header.csv`;

      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "smoothed-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[3]==1)
    {
	    my @args = split(/,/,$argsVec[3]);
    	my $arg0 = $args[0];
    	my $arg1 = $args[1];
        my $now1 = time();
        print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_ab_pssm_no_header.csv -t ab_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";
        $errorMsg = tee_merged {
		`sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_ab_pssm_no_header.csv -t ab_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }

        print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_ab_pssm_no_header.csv -p ab_pssm -n 400 -o originalFeatures/$jobName/$jobName\_ab_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_ab_pssm_no_header.csv -p ab_pssm -n 400 -o originalFeatures/$jobName/$jobName\_ab_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }

      	print "sudo rm originalFeatures/$jobName/$jobName\_ab_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_ab_pssm_no_header.csv`;

        my $now2 = time();
        my $minutesInterval = &getMinutes($now1, $now2);
        print "AB-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[4]==1)
    {
	    my @args = split(/,/,$argsVec[4]);
    	my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
       	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pssm_composition_no_header.csv -t pssm_composition -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pssm_composition_no_header.csv -t pssm_composition -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
       	print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pssm_composition_no_header.csv -p pssm_composition -n 400 -o originalFeatures/$jobName/$jobName\_pssm_composition.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pssm_composition_no_header.csv -p pssm_composition -n 400 -o originalFeatures/$jobName/$jobName\_pssm_composition.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_pssm_composition_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_pssm_composition_no_header.csv`;

       	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
 	 	print "PSSM-composition consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[5]==1)
    {
	    my @args = split(/,/,$argsVec[5]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
       	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_rpm_pssm_no_header.csv -t rpm_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_rpm_pssm_no_header.csv -t rpm_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
       	print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_rpm_pssm_no_header.csv -p rpm_pssm -n 400 -o originalFeatures/$jobName/$jobName\_rpm_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_rpm_pssm_no_header.csv -p rpm_pssm -n 400 -o originalFeatures/$jobName/$jobName\_rpm_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_rpm_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_rpm_pssm_no_header.csv`;

       	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "RPM-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[6]==1)
    {
	    my @args = split(/,/,$argsVec[6]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
      	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_s_fpssm_no_header.csv -t s_fpssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_s_fpssm_no_header.csv -t s_fpssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_s_fpssm_no_header.csv -p s_fpssm -n 400 -o originalFeatures/$jobName/$jobName\_s_fpssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_s_fpssm_no_header.csv -p s_fpssm -n 400 -o originalFeatures/$jobName/$jobName\_s_fpssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_s_fpssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_s_fpssm_no_header.csv`;

      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "S-FPSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[7]==1)
    {
	    my @args = split(/,/,$argsVec[7]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
       	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_dpc_pssm_no_header.csv -t dpc_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_dpc_pssm_no_header.csv -t dpc_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
       	print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_dpc_pssm_no_header.csv -p dpc_pssm -n 400 -o originalFeatures/$jobName/$jobName\_dpc_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_dpc_pssm_no_header.csv -p dpc_pssm -n 400 -o originalFeatures/$jobName/$jobName\_dpc_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_dpc_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_dpc_pssm_no_header.csv`;

       	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "DPC-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[8]==1)
    {
	    my @args = split(/,/,$argsVec[8]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
   		print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm_no_header.csv -t k_separated_bigrams_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm_no_header.csv -t k_separated_bigrams_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm_no_header.csv -p k_separated_bigrams_pssm -n 400 -o originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm_no_header.csv -p k_separated_bigrams_pssm -n 400 -o originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_k_separated_bigrams_pssm_no_header.csv`;

      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "k-separated-bigrams-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[9]==1)
    {
	    my @args = split(/,/,$argsVec[9]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
   		print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_tri_gram_pssm_no_header.csv -t tri_gram_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_tri_gram_pssm_no_header.csv -t tri_gram_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_tri_gram_pssm_no_header.csv -p tri_gram_pssm -n 8000 -o originalFeatures/$jobName/$jobName\_tri_gram_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_tri_gram_pssm_no_header.csv -p tri_gram_pssm -n 8000 -o originalFeatures/$jobName/$jobName\_tri_gram_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_tri_gram_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_tri_gram_pssm_no_header.csv`;

      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "tri-gram-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[10]==1)
    {
	    my @args = split(/,/,$argsVec[10]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
   		print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_eedp_no_header.csv -t eedp -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_eedp_no_header.csv -t eedp -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_eedp_no_header.csv -p eedp -n 400 -o originalFeatures/$jobName/$jobName\_eedp.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_eedp_no_header.csv -p eedp -n 400 -o originalFeatures/$jobName/$jobName\_eedp.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_eedp_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_eedp_no_header.csv`;
      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "EEDP consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[11]==1)
    {
	    my @args = split(/,/,$argsVec[11]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
  	 	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_tpc_no_header.csv -t tpc -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_tpc_no_header.csv -t tpc -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_tpc_no_header.csv -p tpc -n 400 -o originalFeatures/$jobName/$jobName\_tpc.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_tpc_no_header.csv -p tpc -n 400 -o originalFeatures/$jobName/$jobName\_tpc.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_tpc_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_tpc_no_header.csv`;
      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "TPC consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[12]==1)
    {
	    my @args = split(/,/,$argsVec[12]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
   		print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_edp_no_header.csv -t edp -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_edp_no_header.csv -t edp -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_edp_no_header.csv -p edp -n 20 -o originalFeatures/$jobName/$jobName\_edp.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_edp_no_header.csv -p edp -n 20 -o originalFeatures/$jobName/$jobName\_edp.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }

      	print "sudo rm originalFeatures/$jobName/$jobName\_edp_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_edp_no_header.csv`;
   		my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "EDP consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[13]==1)
    {
	    my @args = split(/,/,$argsVec[13]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
   		print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_rpssm_no_header.csv -t rpssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_rpssm_no_header.csv -t rpssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_rpssm_no_header.csv -p rpssm -n 110 -o originalFeatures/$jobName/$jobName\_rpssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_rpssm_no_header.csv -p rpssm -n 110 -o originalFeatures/$jobName/$jobName\_rpssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_rpssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_rpssm_no_header.csv`;
    	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "RPSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[14]==1)
    {
	    my @args = split(/,/,$argsVec[14]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
       	my $now1 = time();
   		print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pse_pssm_no_header.csv -t pse_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pse_pssm_no_header.csv -t pse_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pse_pssm_no_header.csv -p pse_pssm -n 40 -o originalFeatures/$jobName/$jobName\_pse_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pse_pssm_no_header.csv -p pse_pssm -n 40 -o originalFeatures/$jobName/$jobName\_pse_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_pse_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_pse_pssm_no_header.csv`;
      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "Pse-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[15]==1)
    {
	    my @args = split(/,/,$argsVec[15]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
      	my $dimension15 = ($arg0+1)*40;
    	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_dp_pssm_no_header.csv -t dp_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_dp_pssm_no_header.csv -t dp_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_dp_pssm_no_header.csv -p dp_pssm -n $dimension15 -o originalFeatures/$jobName/$jobName\_dp_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_dp_pssm_no_header.csv -p dp_pssm -n $dimension15 -o originalFeatures/$jobName/$jobName\_dp_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_dp_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_dp_pssm_no_header.csv`;
      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "DP-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[16]==1)
    {
	    my @args = split(/,/,$argsVec[16]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
      	my $dimension16 = $arg0*20;
    	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pssm_ac_no_header.csv -t pssm_ac -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pssm_ac_no_header.csv -t pssm_ac -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pssm_ac_no_header.csv -p pssm_ac -n $dimension16 -o originalFeatures/$jobName/$jobName\_pssm_ac.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pssm_ac_no_header.csv -p pssm_ac -n $dimension16 -o originalFeatures/$jobName/$jobName\_pssm_ac.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_pssm_ac_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_pssm_ac_no_header.csv`;
   		my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "PSSM-AC consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[17]==1)
    {
	    my @args = split(/,/,$argsVec[17]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
      	my $dimension17 = $arg0*380;
    	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pssm_cc_no_header.csv -t pssm_cc -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_pssm_cc_no_header.csv -t pssm_cc -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pssm_cc_no_header.csv -p pssm_cc -n $dimension17 -o originalFeatures/$jobName/$jobName\_pssm_cc.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_pssm_cc_no_header.csv -p pssm_cc -n $dimension17 -o originalFeatures/$jobName/$jobName\_pssm_cc.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_pssm_cc_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_pssm_cc_no_header.csv`;
      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "PSSM-CC consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[18]==1)
    {
	    my @args = split(/,/,$argsVec[18]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
    	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_aadp_pssm_no_header.csv -t aadp_pssm -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_aadp_pssm_no_header.csv -t aadp_pssm -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_aadp_pssm_no_header.csv -p aadp_pssm -n 420 -o originalFeatures/$jobName/$jobName\_aadp_pssm.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_aadp_pssm_no_header.csv -p aadp_pssm -n 420 -o originalFeatures/$jobName/$jobName\_aadp_pssm.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo rm originalFeatures/$jobName/$jobName\_aadp_pssm_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_aadp_pssm_no_header.csv`;
   		my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "AADP-PSSM consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[19]==1)
    {
	    my @args = split(/,/,$argsVec[19]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
    	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_aatp_no_header.csv -t aatp -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_aatp_no_header.csv -t aatp -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
   		print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_aatp_no_header.csv -p aatp -n 420 -o originalFeatures/$jobName/$jobName\_aatp.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_aatp_no_header.csv -p aatp -n 420 -o originalFeatures/$jobName/$jobName\_aatp.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
  		print "sudo rm originalFeatures/$jobName/$jobName\_aatp_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_aatp_no_header.csv`;
   		my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "AATP consuming time is $minutesInterval minutes \n";
    }
    if($featureVec[20]==1)
    {
	    my @args = split(/,/,$argsVec[20]);
	    my $arg0 = $args[0];
	    my $arg1 = $args[1];
      	my $now1 = time();
      	print "sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_medp_no_header.csv -t medp -p pssm/$jobName -a $arg0 -b $arg1\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/possum.py -i sequence/$jobName.fasta -o originalFeatures/$jobName/$jobName\_medp_no_header.csv -t medp -p pssm/$jobName -a $arg0 -b $arg1`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
      	print "sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_medp_no_header.csv -p medp -n 420 -o originalFeatures/$jobName/$jobName\_medp.csv\n";

        $errorMsg = tee_merged {
            `sudo python featureCalcCode/headerHandler.py -i originalFeatures/$jobName/$jobName\_medp_no_header.csv -p medp -n 420 -o originalFeatures/$jobName/$jobName\_medp.csv`;
        };

        if( ( $? >> 8 ) != 0 ) {
            $retStr = $retStr."error;";
            $retStr = $retStr.$jobName;
            $retStr = $retStr.";";
            $retStr = $retStr.$selectedFeatures;
            $retStr = $retStr.";";
            $retStr = $retStr.$errorMsg;
            return $retStr;
        }
  		print "sudo rm originalFeatures/$jobName/$jobName\_medp_no_header.csv\n";
  		`sudo rm originalFeatures/$jobName/$jobName\_medp_no_header.csv`;

      	my $now2 = time();
      	my $minutesInterval = &getMinutes($now1, $now2);
      	print "MEDP consuming time is $minutesInterval minutes \n";
    }

    print "Success to generate original features\n";
    my $finishCalculatingTime = time();

    my $calculateInterval = &getMinutes($beginCalculatingTime, $finishCalculatingTime);
    print "Calculating features consuming time is $calculateInterval minutes\n";
    # reorganize result files
    if(index($retStr,"warning") != -1)
    {
        io("$tomcat_storage_dir/$jobName/$jobName\_pssm_missed_sequences.txt")->print($sequences_missed_pssm);
    }

    my $beginZipTime = time();
    print "sudo zip $tomcat_storage_dir/$jobName/$jobName\_pssm_files\.zip pssm/$jobName -q -r\n";

    $errorMsg = tee_merged {
        `sudo zip $tomcat_storage_dir/$jobName/$jobName\_pssm_files\.zip pssm/$jobName -q -r`;
    };

    if( ( $? >> 8 ) != 0 ) {
        $retStr = $retStr."error;";
        $retStr = $retStr.$jobName;
        $retStr = $retStr.";";
        $retStr = $retStr.$selectedFeatures;
        $retStr = $retStr.";";
        $retStr = $retStr.$errorMsg;
        return $retStr;
    }
    print "sudo zip $tomcat_storage_dir/$jobName/$jobName\_pssm_features\.zip originalFeatures/$jobName -q -r\n";

    $errorMsg = tee_merged {
        `sudo zip $tomcat_storage_dir/$jobName/$jobName\_pssm_features\.zip originalFeatures/$jobName -q -r`;
    };

    if( ( $? >> 8 ) != 0 ) {
        $retStr = $retStr."error;";
        $retStr = $retStr.$jobName;
        $retStr = $retStr.";";
        $retStr = $retStr.$selectedFeatures;
        $retStr = $retStr.";";
        $retStr = $retStr.$errorMsg;
        return $retStr;
    }
    print "sudo mv originalFeatures/$jobName $tomcat_storage_dir/$jobName/$jobName\_features\n";

    $errorMsg = tee_merged {
        `sudo mv originalFeatures/$jobName $tomcat_storage_dir/$jobName/$jobName\_features`;
    };

    if( ( $? >> 8 ) != 0 ) {
        $retStr = $retStr."error;";
        $retStr = $retStr.$jobName;
        $retStr = $retStr.";";
        $retStr = $retStr.$selectedFeatures;
        $retStr = $retStr.";";
        $retStr = $retStr.$errorMsg;
        return $retStr;
    }
    my $finishZipTime = time();
    my $zipInterval = &getMinutes($beginZipTime, $finishZipTime);
    print "Compressing features consuming time is $zipInterval minutes\n";

    my $finishWorkflowTime = time();
    my $totalInterval = &getMinutes($beginWorkflowTime, $finishWorkflowTime);
    print "Total consuming time is $totalInterval minutes\n";
    print "Finish the workflow\n";
    #########################################################
    ################### Need modify   #######################
    #########################################################


    # integrate results
    if(index($retStr,"warning") == -1)
    {
        $retStr = $retStr."results;";
    }
    $retStr = $retStr.$jobName;
    $retStr = $retStr.";";
    $retStr = $retStr.$selectedFeatures;
    $retStr = $retStr.";";
    if(index($retStr,"warning") == -1)
    {
        $retStr = $retStr."No results";
    }
    elsif(index($retStr,"warning") != -1){
        $retStr = $retStr.$missed_pssm_count;
    }

    #my @lines = io("predictResults/$jobName\_results.csv")->slurp;
    #print $retStr;
    #io("$jobName\_ret.txt")->print($retStr);

    ##################################################
    ####     Begin to clean intermediate files    ####
    ##################################################
    print "Begin to clean intermediate files: \n";
    #print "sudo rm -f sequence/$jobName.fasta\n";
    #`sudo rm -f sequence/$jobName.fasta`;
    print "sudo rm -fr seperatedSequence/$jobName\n";
    `sudo rm -fr seperatedSequence/$jobName`;
    print "sudo rm -fr pssm/$jobName\n";
    `sudo rm -fr pssm/$jobName`;

    print "Success to clean intermediate files!\n";
    ############################################################
    return $retStr;
}

sub getTime
{
    my $time = shift || time();
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);

    $year += 1900;
    $mon ++;

    $min  = '0'.$min  if length($min)  < 2;
    $sec  = '0'.$sec  if length($sec)  < 2;
    $mon  = '0'.$mon  if length($mon)  < 2;
    $mday = '0'.$mday if length($mday) < 2;
    $hour = '0'.$hour if length($hour) < 2;

    my $weekday = ('Sun','Mon','Tue','Wed','Thu','Fri','Sat')[$wday];

    return { 'second' => $sec,
             'minute' => $min,
             'hour'   => $hour,
             'day'    => $mday,
             'month'  => $mon,
             'year'   => $year,
             'weekNo' => $wday,
             'wday'   => $weekday,
             'yday'   => $yday,
             'date'   => "$year-$mon-$mday"
          };
}

sub getMinutes
{
  my ($beginTime, $finishTime) = @_;
  my $timeInterval = $finishTime - $beginTime;
  my $minutesInterval = $timeInterval/60.0;
  my $formatInterval = sprintf("%.3f", $minutesInterval);
  return $formatInterval;
}

sub getSeqByPath{
    my ($seq_path) = @_;
    my $objs = Bio::SeqIO->new(-file=>$seq_path, -format=>'fasta');
    my $obj = $objs->next_seq;
    return $obj->seq;
}
sub getTitleAndSeqByPath{
    my ($seq_path) = @_;
    my $objs = Bio::SeqIO->new(-file=>$seq_path, -format=>'fasta');
    my $seq_obj = $objs->next_seq;
    my $display_name = $seq_obj->display_name;
    my $desc = $seq_obj->desc;
    my $seq = $seq_obj->seq;
    my $seq_type = $seq_obj->alphabet;
    my $seq_length = $seq_obj->length;
    my $seqContent = ">".$display_name." ".$desc."\n".$seq."\n";
    return $seqContent;
}

sub updateStatus {
    print "start to update status in job table!\n";

    my ($jobName,$status,$timeType)=@_;
    my $timeHandler= $timeType==1?"startTime":"endTime";
    print '$jobName =';
    print "$jobName\n";
    print '$timeHandler =';
    print "$timeHandler\n";
    my $dbh = DBI->connect("DBI:mysql:possum:127.0.0.1", 'root', 'POSSUMdb');

    my $sql=qq{UPDATE job SET   status = $status, $timeHandler = now()  WHERE jobName = '$jobName'};
    print "sql statements to be execute:\n";
    print "$sql"."\n";

    my $sth = $dbh->prepare($sql);
    $sth->execute() or die $DBI::errstr;

    $sth->finish();
    $dbh->commit;
    $dbh->disconnect();
    print "Success to update status in job table!\n";
}

sub getCacheFileName{
    my ($table,$seq)=@_;
    print '$table =';
    print "$table\n";
    print '$seq =';
    print "$seq\n";
    my $dbh = DBI->connect("DBI:mysql:possum:127.0.0.1", 'root', 'POSSUMdb');
    my $sql =  qq{SELECT output FROM $table WHERE seq = '$seq'};
    print "sql statements to be execute:\n";
    print "$sql"."\n";
    my $sth = $dbh->prepare($sql);
    $sth->execute() or die $DBI::errstr;
    my @row = $sth->fetchrow_array();

    $sth->finish();
    $dbh->commit;
    $dbh->disconnect();

    return $row[0];
}

sub storeCacheFileName{
    my ($table,$seq,$output)=@_;
    print '$table =';
    print "$table\n";
    print '$seq =';
    print "$seq\n";
    print '$output =';
    print "$output\n";
    my $dbh = DBI->connect("DBI:mysql:possum:127.0.0.1", 'root', 'POSSUMdb');
    my $sql =  qq{insert into $table (seq,output) values('$seq','$output')};
    print "sql statements to be execute:\n";
    print "$sql"."\n";
    my $sth = $dbh->prepare($sql);
    $sth->execute() or die $DBI::errstr;

    $sth->finish();
    $dbh->commit;
    $dbh->disconnect();
}
