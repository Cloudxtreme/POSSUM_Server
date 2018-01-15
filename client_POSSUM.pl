#!/usr/bin/perl
use SOAP::Transport::HTTP;
SOAP::Transport::HTTP::CGI
-> dispatch_to('POSSUM')
-> handle;


package POSSUM;

#use strict;
use warnings;
use Gearman::Client;
use Storable;
use Storable qw(freeze);
use Storable qw(thaw);
use DBI;
#use MIME::Lite;
#use Mail::Sendmail;
use Email::MIME;
use Email::Sender::Simple qw(sendmail);
use Email::Sender::Transport::SMTP qw();
use Try::Tiny;

sub testPredict
{
	return "Success!";
}

sub predict
{
	# The first parameter class is a default value of this file, you just need to receive it without solving it.
	local ($class,$jobName,$seqs,$seqNum,$selectedFeatures,$email,$databaseType,$arguments,$iterations,$E_value)=@_;
	local $retStr="";
	$seqs =~ s/&gt;/>/g;
	#return $seqs;

	# fork this process
	my $pid = fork();
	if ($pid == 0)
	{
		# do this in the child, or apache server will response to the web request after child process is finished,
		# which will lead to timeout when the child process is high time-consuming.
		open STDIN, "</dev/null";
		open STDOUT, ">/dev/null";
		open STDERR, ">/dev/null";

		print "start new client \n";
		my $client = Gearman::Client->new;
		print "finish new client \n";

		print "start job_servers \n";
		$client->job_servers('127.0.0.1',4730);
		print "finish job_servers \n";

		print "start new_task_set \n";
		my $tasks = $client->new_task_set;
		print "finish new_task_set \n";

		print "start add_task \n";

		my @rows=($class,$jobName,$seqs,$seqNum,$selectedFeatures,$email,$databaseType,$arguments,$iterations,$E_value);
		$tasks->add_task(
		      generatePSSM => freeze(\@rows),

		      { on_complete => \&complete },
		      );
		print "finish add_task \n";

		print "client start wait \n";
		$tasks->wait;
		print "client finish wait \n";
		exit;
	}
	$retStr = $retStr."Request has been accepted and is running background\n";
	return $retStr;
}


sub complete{

    my $ret = ${ $_[0] };
    #io("complete.txt")->print($ret);
    print $ret, "\n";
    #######################################updatePredictResults##############
    &updatePredictResults($ret);

	my @rets = split(/;/,$ret);
	my $retMarker = $rets[0];
	my $jobResults = $rets[3];

    if($email ne "")
    {
		if($retMarker eq "error") {
			&updatePredictResultsTable($jobName,$jobResults);
			&POSSUM_send_error_mail($jobName,$email,$jobResults);
		} else {
			&POSSUM_send_mail($jobName,$email);
		}

    }

}

sub updatePredictResultsTable {
	print "start to update status in predictresult table!\n";

	my ($jobName,$predictionResult)=@_;

	print '$jobName =';
	print "$jobName\n";
	print '$predictionResult =';
	print "$predictionResult\n";
	my $dbh = DBI->connect("DBI:mysql:possum:127.0.0.1", 'root', 'POSSUMdb');
	my $sql=qq{UPDATE computeresult SET   resultInfo = "$predictionResult"  WHERE jobName = '$jobName'};
	print "sql statements to be execute:\n";
	print "$sql"."\n";

	my $sth = $dbh->prepare($sql);
	$sth->execute() or die $DBI::errstr;

	$sth->finish();
	$dbh->commit;
	$dbh->disconnect();
	print "Success to update status in predictresult table!\n";
}

sub updatePredictResults {
	print "start to insert record in predictresult table!\n";

	my ($ret) = @_;
	print '$ret =';
	print "$ret\n";
	my $STATUS_COMPLETE = 2;
	my $STATUS_ERROR = 3;
	my $STATUS_WARNING = 4;
	my $END_TIME = 2;
	my @rets = split(/;/,$ret);
	my $retMarker = $rets[0];
	my $jobName = $rets[1];
	my $jobResults = $rets[3];
	if($retMarker eq "results")
	{
	    #update status to 2
	    &updateStatus($jobName, $STATUS_COMPLETE, $END_TIME,$jobResults);
	}
	elsif($retMarker eq "warning"){
		#update status to 4
    	&updateStatus($jobName, $STATUS_WARNING, $END_TIME,$jobResults);
	}
	else
	{
	    #update status to 3
	    &updateStatus($jobName, $STATUS_ERROR, $END_TIME,$jobResults);
	    #update predictResults;
      	    #&updatePredictResultsTable($jobName,$jobResults);a
	}


	print "Success to insert record in predictresult table!\n";
}

sub updateStatus {
	print "start to update status in job table!\n";

	my ($jobName,$status,$timeType,$jobResults)=@_;
	my $timeHandler= $timeType==1?"startTime":"endTime";
	my $pssmMissedCount = ($status==4?$jobResults:0);
	print '$jobName =';
	print "$jobName\n";
	print '$timeHandler =';
	print "$timeHandler\n";
	my $dbh = DBI->connect("DBI:mysql:possum:127.0.0.1", 'root', 'POSSUMdb');

	my $sql=qq{UPDATE job SET   status = $status, $timeHandler = now(),pssmMissedCount = $pssmMissedCount WHERE jobName = '$jobName'};
	print "sql statements to be execute:\n";
	print "$sql"."\n";

	my $sth = $dbh->prepare($sql);
	$sth->execute() or die $DBI::errstr;

	$sth->finish();
	$dbh->commit;
	$dbh->disconnect();
	print "Success to update status in job table!\n";
}

sub POSSUM_send_mail
{
    my ($jobName,$to) = @_;

    print "Email Sent to $to\n";
    my $subject = "POSSUM: Your job $jobName is finished!";
    my $mailContent = "Dear Sir/Madam:\n\n";
    $mailContent = $mailContent."Thanks for using POSSUM. ";
    $mailContent = $mailContent."Your job $jobName is finished!\n";
    $mailContent = $mailContent."Please use the following web link to access the results:\n";
    $mailContent = $mailContent."http://possum.erc.monash.edu/showJobDetail.action?jobName=$jobName";
    my $message = Email::MIME->create(
        header_str => [
            From    => 'admin@possum.erc.monash.edu.au',
            To      => $to,
            Subject => $subject,
        ],
        attributes => {
            encoding => 'quoted-printable',
            charset  => "US-ASCII",
        },
        body_str => $mailContent,
    );


    try {
        sendmail(
            $message,
            {
                From      => "admin\@possum.erc.monash.edu.au",
                To        => $to,
                transport => Email::Sender::Transport::SMTP->new({
                    host  => "smtp.monash.edu",
                    port  => 25,
                })
            }
        );
    } catch {
        warn "sending failed: $_";
    };
    print "Email Sent Successfully\n";
}

sub POSSUM_send_error_mail
{
	my ($jobName,$to,$errorInfo) = @_;

	$subject = "POSSUM: Your job $jobName went wrong!";
	my $mailContent = "Dear Sir/Madam:\n\n";
	$mailContent = $mailContent."Thanks for using POSSUM. \n";
	$mailContent = $mailContent."An unexpected error occurs during the process of feature generating, please refer to the following information for details:\n";
	$mailContent = $mailContent.$errorInfo;
	$mailContent = $mailContent."\n";
	$mailContent = $mailContent."This error can be caused by several possible reasons:\n";
	$mailContent = $mailContent."1. Your input file was incorrectly formatted, which failed to meet the requirment of the POSSUM web server.\n";
	$mailContent = $mailContent."2. The POSSUM web server crashed due to the unknown system error.\n";
	$mailContent = $mailContent."\n";
	$mailContent = $mailContent."To tackle this issue, you'd better\n";
	$mailContent = $mailContent."1. Format your input file in accordance with the input formats described in the help page http://possum.erc.monash.edu/help.jsp#input and try again http://possum.erc.monash.edu/server.jsp.\n";
	$mailContent = $mailContent."2. Send an email to Chris (chris\@nohup.cc) or Young (young\@nohup.cc) together with your input file and job url http://possum.erc.monash.edu/showJobDetail.action?jobName=$jobName, if this issue still can't be tackled.\n";
	$mailContent = $mailContent."\n";

	my $message = Email::MIME->create(
        header_str => [
            From    => 'admin@possum.erc.monash.edu.au',
            To      => $to,
            Subject => $subject,
        ],
        attributes => {
            encoding => 'quoted-printable',
            charset  => "US-ASCII",
        },
        body_str => $mailContent,
    );


    try {
        sendmail(
            $message,
            {
                From      => "admin\@possum.erc.monash.edu.au",
                To        => $to,
                transport => Email::Sender::Transport::SMTP->new({
                    host  => "smtp.monash.edu",
                    port  => 25,
                })
            }
        );
    } catch {
        warn "sending failed: $_";
    };
    print "Email Sent Successfully\n";
}

sub POSSUM_send_mail_used_before
{
	my ($jobName,$to) = @_;

	$subject = "POSSUM: Your job $jobName is finished!";
	my $mailContent = "Dear Sir/Madam:\n\n";
	$mailContent = $mailContent."Thanks for using POSSUM. ";
	$mailContent = $mailContent."Your job $jobName is finished!\n";
	$mailContent = $mailContent."Please use the following web link to access the results:\n";
	$mailContent = $mailContent."http://possum.erc.monash.edu/showJobDetail.action?jobName=$jobName";

	%mail = ( To => $to,From => 'admin@possum.erc.monash.edu',Subject => $subject, Message => $mailContent );
	sendmail(%mail) or die $Mail::Sendmail::error;
	print "Email Sent Successfully\n";
}

sub POSSUM_send_error_mail_used_before
{
	my ($jobName,$to,$errorInfo) = @_;

	$subject = "POSSUM: Your job $jobName went wrong!";
	my $mailContent = "Dear Sir/Madam:\n\n";
	$mailContent = $mailContent."Thanks for using POSSUM. \n";
	$mailContent = $mailContent."An unexpected error occurs during the process of feature generating, please refer to the following information for details:\n";
	$mailContent = $mailContent.$errorInfo;
	$mailContent = $mailContent."\n";
	$mailContent = $mailContent."This error can be caused by several possible reasons:\n";
	$mailContent = $mailContent."1. Your input file was incorrectly formatted, which failed to meet the requirment of the POSSUM web server.\n";
	$mailContent = $mailContent."2. The POSSUM web server crashed due to the unknown system error.\n";
	$mailContent = $mailContent."\n";
	$mailContent = $mailContent."To tackle this issue, you'd better\n";
	$mailContent = $mailContent."1. Format your input file in accordance with the input formats described in the help page http://possum.erc.monash.edu/help.jsp#input and try again http://possum.erc.monash.edu/server.jsp.\n";
	$mailContent = $mailContent."2. Send an email to Chris (chris\@nohup.cc) or Young (young\@nohup.cc) together with your input file and job url http://possum.erc.monash.edu/showJobDetail.action?jobName=$jobName, if this issue still can't be tackled.\n";
	$mailContent = $mailContent."\n";

	# $mailContent = $mailContent."Please use the following web link to access the results:\n";
	# $mailContent = $mailContent."http://possum.erc.monash.edu/showJobDetail.action?jobName=$jobName";

	%mail = ( To => $to,From => 'admin@possum.erc.monash.edu',Subject => $subject, Message => $mailContent );
	sendmail(%mail) or die $Mail::Sendmail::error;
	print "Email Sent Successfully\n";
}
