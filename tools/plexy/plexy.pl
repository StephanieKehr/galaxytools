#!/usr/bin/env perl

#Stephanie Kehr

#my $PATH = "/homes/bierdepot/steffi/bin/";
##my $PATH = ".";
#if(!-e $PATH."/RNAplex"){
#    print STDERR "RNAplex is not installed in the given PATH\n";exit(0);
#}

use warnings;
use strict;
use Getopt::Std;


my %opts;
getopts('f:t:T:p:o:e:l',\%opts);

if(defined($opts{f}) && defined($opts{o}) && (defined($opts{t}) || defined($opts{T})) ) {
    
    #default energy-treshold
    if(!defined($opts{e})){ $opts{e} = -7.70};
    my $targetLoc;
    if(!defined($opts{t})){ $targetLoc = $opts{T}; } else { $targetLoc = $opts{t}; }
    
    #read fa-file
    my $fasta = &readFa($opts{f});

    #user-input
    my $outDir = $opts{o};
    $outDir .= "/" if($outDir !~ /\/$/);
    #my $targetDir = $opts{t};

    #create output-dir
    if(!-e $outDir){
	`mkdir $outDir`;
    }

    foreach my $sno(@{$fasta}){    
	my ($boxStart, $boxPrimeStart) ;
	my ($box, $boxPrime);
	my ($snoID,$targetID);
	my ($head,$seq);

	my @tmp = split(/\n/,$sno);
	$head = $tmp[0];
	$seq = $tmp[1];

	@tmp = split(/_/,$head);

	#generalInfo
	($snoID = $tmp[0]."_".$tmp[1]) =~ s/>//;

	#D-box-information
	$boxStart = $tmp[-1];
	$boxStart =~ s/\D+$//;
	$box = $tmp[-2];

	#check whether D-box annotation is correct,stop otherwise
	if(!defined($boxStart) || !defined($box) || $boxStart =~ /\D/ || $box =~ /[^AGCTU]{4}/){
	    print STDOUT $head,"\n",$seq,"\n";
	    print STDERR "WRONG HEADER! CAN'T IDENTIFY D_BOX!\n";
	    exit(0);
	}
	else{
	    #print $head,"\n",$seq,"\n";
	    #compute D-box targets
	    print STDOUT "#D-box targets:\n";

	    &computeTargetCD($snoID,$seq,$boxStart,$targetLoc,$outDir,"D");	


	    #D'-box-information
	    $boxPrimeStart = $tmp[-5];
	    $boxPrime = $tmp[-6];

	    #check whether D'-box annotation is also available
	    if(!defined($boxPrimeStart) || !defined($boxPrime) || $boxPrimeStart =~ /\D/ || $boxPrime =~ /[^AGCTU]{4}/){
		$tmp[0] = "CAN'T IDENTIFY D'_BOX!\n";
	    }
	    else{
	        #compute D'-box targets
	    	print STDOUT "#D'-box targets:\n";

		&computeTargetCD($snoID,$seq,$boxPrimeStart,$targetLoc,$outDir,"D'");
	    }
	}
	
    }
}
else{
    print "USAGE:\n
plexy.pl -o [directory] -f [fa-file] -t|T [directory|File]

\n \t -f : input fasta-file with CD-box snoRNA sequences
\n \t -o : directory for output-file
\n \t -t : directory with putative target RNA sequences (or use -T)
\n \t -T : single file with putative target RNA sequences (or use -t)
\n\n OPTIONAL
\n \t [-l] : option for large datasets. Results are not kept in the memory during runtime
                but are returned directly and the output is reduced. This makes it possible
                to search transcript wide, but the results need some postprocessing as inter-
                actions are not unique and not sorted for mfe values.
\n \t [-p] : directory containing accessibility-profiles
\n \t [-e] : MFE-treshold, only targets with better MFE are shown, default value is -7.70
";

}

########################################################################
sub computeTargetCD
#finds regions which can bind to the InteractionRegion using RNAplex
########################################################################
{
    my ($id,$seq,$Dbox,$targetLoc,$dir,$boxType) = @_;
   
    my $target;
    my @results;
    
    my $ia = &cutInteractionRegion($Dbox,$seq);
    
    #collects all possible targetRNAs
    if($opts{t}){
	my @targetRNAs = glob($targetLoc."/*.fa*");
	
	
	#for each possible targetRNA
	foreach my $tarRNAfile(@targetRNAs){
	
	    #read targetSequence
	    open(TAR,"<".$tarRNAfile);
	    my @tar = <TAR>;
	    close(TAR);
	    
	    #only name of targetRNA, so throw away the path
	    my $tarRNA = $tar[0];
	    if($tarRNA !~ /^>/){
		print STDERR "File containing putative target sequence: ",$tarRNAfile," is not fasta-formatted!\nContinue with next putative target\n";
		next;
	    }
	    my $tarSeq = $tar[1];
	    chomp($tarSeq);
	    
	    #print file for RNAplex
	    open(PLEX, ">".$dir."$$\_plex.in");
	    print PLEX $tarRNA,$tarSeq,"\n>InteractionRegion\n",$ia,"\n";
	    close(PLEX);
	    
	    
	    chomp($tarRNA);
	    $tarRNA =~ s/>//;
	    
	    #$tarRNAfile = $tarRNA;
	    
	    my @hits = '';
	    if(defined($opts{p})){
		
		#no dummy-profile
		if(defined($opts{p}) && !-e $opts{p}."InteractionRegion_openen"){
		    &writeDummyProfile;
		}
		
		#no profile for targetRNA
		if(!-e $opts{p}."/".$tarRNA."_openen"){
		    print STDERR "NO PROFILE FOUND FOR ",$tarRNAfile,".fa IN ",$opts{p},"\nPLEASE RUN RNAplfold FROM THE VIENNA RNA PACKAGE\n";
		    exit(0); 
		} 
		
		#@hits = `$PATH/RNAplex <$dir$$\_plex.in -l 19 -e $opts{e} -a $opts{p} -f 2`;
		@hits = `RNAplex <$dir$$\_plex.in -l 19 -e $opts{e} -a $opts{p} -f 2`;
	    }
	    else{
		#@hits = `$PATH/RNAplex <$dir$$\_plex.in -l 20 -e $opts{e} -f 2`;
		@hits = `RNAplex <$dir$$\_plex.in -l 20 -e $opts{e} -f 2`;
	    }
	    
	    push(@results, @{&filterHits(\@hits,$tarRNA,$tarSeq,$ia,$id,$boxType)});
	}
    }
    #single targetFile is defined with opts{T}
    else{
	open(TAR,"<".$targetLoc);
	my @tar = <TAR>;
	close(TAR);

	for(my $i = 0 ; $i < scalar(@tar) ; $i+=2){
	    my $tarRNA = $tar[$i];

	    if($tarRNA !~ /^>/){
		print STDERR "File containing putative target sequences: ",$targetLoc," is not fasta-formatted!\n";
		exit(0);
	    }
	    
	    my $tarSeq = $tar[$i+1];
	    chomp($tarSeq);
	    
	    #print file for RNAplex
	    open(PLEX, ">".$dir."$$\_plex.in");
	    print PLEX $tarRNA,$tarSeq,"\n>InteractionRegion\n",$ia,"\n";
	    close(PLEX);

	    chomp($tarRNA);
	    $tarRNA =~ s/>//; 

	    my @hits = '';
	    if(defined($opts{p})){
		
		#no dummy-profile
		if(defined($opts{p}) && !-e $opts{p}."InteractionRegion_openen"){
		    &writeDummyProfile;
		}
		
		#no profile for targetRNA
		if(!-e $opts{p}."/".$tarRNA."_openen"){
		    print STDERR "NO PROFILE FOUND FOR ",$tarRNA,".fa IN ",$opts{p},"\nPLEASE RUN RNAplfold FROM THE VIENNA RNA PACKAGE\n";
		    exit(0); 
		} 
		
		#@hits = `$PATH/RNAplex <$dir$$\_plex.in -l 19 -e $opts{e} -a $opts{p} -f 2`;
		@hits = `RNAplex <$dir$$\_plex.in -l 19 -e $opts{e} -a $opts{p} -f 2`;
	    }
	    else{
		#@hits = `$PATH/RNAplex <$dir$$\_plex.in -l 20 -e $opts{e} -f 2`;
		@hits = `RNAplex <$dir$$\_plex.in -l 20 -e $opts{e} -f 2`;
	    }

	    push(@results, @{&filterHits(\@hits,$tarRNA,$tarSeq,$ia,$id,$boxType)});

	}
    }
    
    #sort and unique in case option l is not defined
    if(scalar(@results)>0){

    	#sort results after energy-values
    	my @sorted = sort {$a->{energy} <=> $b->{energy}} @results;
    	@results = '';

    	#make modificated sites unique
    	my @uniq;
    	my %seen = ();
    	foreach (@sorted) {
		my $tmp = $_;
		#there are some variants for the snRNAs, but only the best for each snRNA should be regarded
		(my $tmptar = $tmp->{targetRNA}) =~ s/\.\d+$//;
		my $item = $tmptar."-".$tmp->{mod};
	
		push(@uniq,$tmp) unless $seen{$item}++;
    	}
    	@sorted = '';
    
    	my $targets;
    	if(scalar(@uniq) > 0 ){
		foreach(@uniq){
	    	$target = $_;
		    printf STDOUT "%-10s\t%-15s\t%-20s\t%-45s%-45s\n",$id."\t".$boxType,$_->{targetRNA}."-".$_->{mod},$_->{energy},$_->{structure},$_->{seq};
		}
    	}
    	else{
	    printf STDOUT "No targets found!\n";
    	}
    }
}

########################################################################
sub hasNoBulge
#checks weather the duplex has an bulge
#returns 1 if there is no bulge and -1 otherwise
########################################################################
{
    my($struc1,$struc2) = @_;
    my $bool = 1;

    my $offset = 0;
    my $mismatch = index($struc1, ".", $offset);
    #for every occurence of a mismatch(dot)
    while ($mismatch != -1) {
	#there has to be a mismatch in the other sequence,too			
	if(substr(reverse($struc2),$mismatch,1) ne "."){
	    last;
	}
	$offset = $mismatch + 1;
	$mismatch = index($struc1, ".", $offset);
    }
    #jumped out of while-loop by last (otherwise mismatch would be -1), so a bulge was found
    if($mismatch != -1){
	$bool = -1;
    }
    
    return $bool;
}

#######################################################################
sub filterHits
#filters the duplex gained by RNAplex for additonal rules
#######################################################################
{
    my ($hits,$tarRNA,$tarSeq,$ia,$id,$boxType) = @_;
    my @results;
 
   foreach my $hit(@{$hits}){
	chomp($hit);
	#print $hit,"\n";
	
	#not needed output of RNAplex
	if($hit =~ /^>/ || $hit =~ /^$/){
	    next;
	}
	#the duplexes
	else{
	    
	    my $hash;
	    my @parts = split(/\s+/,$hit);
	    $hash->{targetRNA} = $tarRNA;
	    if($parts[4] =~ /\-/){
		($hash->{energy} = $parts[4]) =~ s/[\(\)]//g;
	    }
	    else{
		$hash->{energy} = $parts[5] =~ s/[\(\)]//g;
	    }
	    #recomputed energy value from RNAplex
	    if(defined $opts{p}){
	    	$hit =~ /^[^\[]+\[(\-\d+.\d+)/;
	    	$hash->{energy} = $1;
	    }
	    $hash->{targetPos} = $parts[1];
	    $hash->{bindPos} = $parts[3]; 
	    $hash->{raw_structure} = $parts[0];
	    #structure of Binding
	    my @struc = split(/&/,$hash->{raw_structure});
	    my $tarStart = (split(/,/,$hash->{targetPos}))[0];
	    my ($duplexStart,$duplexEnd) = ((split(/,/,$hash->{bindPos}))[0],(split(/,/,$hash->{bindPos}))[1]);
	    #RNAplex needs this when used without profiles?
	    if($struc[0] =~ /^\./ && $struc[1] =~ /[^\.]$/){
		$struc[0] =~ s/^\.//;
		$tarStart++;
	    }
	    $hash->{structure} = $struc[0]."&".$struc[1];
	    
	    #duplex has to start at least 2 nts before methylated site
	    if($duplexEnd <= 17){
		#print STDERR "DUPLEX TOO FAR FROM BOX:\t ",$duplexEnd,"\n";
		next;
	    }
	    $hash->{bindSeq} = substr($ia,$duplexStart-1,$duplexEnd - $duplexStart+1);
	    $hash->{tarSeq} = substr($tarSeq,$tarStart-1,(split(/,/,$hash->{targetPos}))[1] - $tarStart+1);
	    my $distToBox = 20 - $duplexEnd;
	    
	    $hash->{mod} = $tarStart + 4 - $distToBox;
	    
	    $hash->{seq} = $hash->{tarSeq}."&".$hash->{bindSeq};
	    
	    #duplex-length has to be at least 7 basepairs long
	    if(length($hash->{tarSeq}) < 7){
		#print STDERR "DUPLEX IS TO SHORT:\t ",length($hash->{tarSeq}),"\n";
		next;
	    }
	    
	    #at beginning and end of the targetSeq a bulge is allowed to occur
	    my $tarStruc = substr($struc[0],1,length($struc[0])-2);
	    my $bindStruc =  substr($struc[1],1,length($struc[1])-2);

	    #throw away duplexes with right or left bulge
	    if(&hasNoBulge($tarStruc,$bindStruc)== -1 || &hasNoBulge($struc[1],$struc[0])== -1){
		#print STDERR "FOUND BULGE IN:\t", $hash->{structure}."\n";
		next;
	    }
	    
	    #check if 5th position upstream of D/D'-box, so methylated nucleotide forms a Watson-Crick basepair
	    $hash->{modnt} = (split(//,$hash->{tarSeq}))[4 - $distToBox];
	    my $ntPair = (split(//,$hash->{bindSeq}))[-5 + $distToBox];
	    if( defined($ntPair) && (($hash->{modnt} eq "A" && $ntPair ne "U" &&  $ntPair ne "T") || (($hash->{modnt} eq "U" || $hash->{modnt} eq "T") && $ntPair ne "A") || ($hash->{modnt} eq "G" && $ntPair ne "C") || ($hash->{modnt} eq "C" && $ntPair ne "G")) ){
		#print STDERR "FOUND NO WATSON-CRICK PAIR:\t", $hash->{modnt}, "--", $ntPair,"\n";
		next;
	    }
	    $hash->{tarSeq} = substr($hash->{tarSeq},0, 5 - $distToBox)."m".substr($hash->{tarSeq},5 - $distToBox, length($hash->{tarSeq}));
	    $hash->{seq} = $hash->{tarSeq}."&".$hash->{bindSeq};
	    
	    #allow only duplexes with at most 1(maybe change to two) mismatches in duplex core-region between 3rd and 11th nt upsteam ot D/D'-box
	    if(grep( /\./,split(//,substr($struc[1],2,9))) > 1){
		#print STDERR "MORE THAN ONE MISMATCH IN DUPLEX CORE REGION:\t",substr($struc[1],2,9) ,"\n";
		next;
	    }
	    #filter hits with energy lower than treshold (needed because of recomputing of RNAplex energys)
	    if($hash->{energy} > $opts{e}){
		#print STDERR "MFE:",$hash->{energy}," IS TOO HIGH\n";
		next;
	    }
	    if(defined($opts{l})){
		print STDOUT $id,"|",$boxType,"|",,$hash->{targetRNA}."-".$hash->{mod},"|",$hash->{energy},"|",$hash->{structure},"|",$hash->{seq},"\n";
	    }	
	    else{
		push(@results,$hash);
	    }	
	}
    }
    return \@results;
}



########################################################################
sub cutInteractionRegion
#cuts the region upstream of the D or D' box, which does potentialy bind
#to RNA region which will be methylated
########################################################################
{
    my ($boxPos,$seq) = @_;
    #cuts 20 nucleotides before D|D'-box from sequence
    my $interaction;
    if($boxPos >  20){
	$interaction = substr($seq,$boxPos-21,20);
    } else{
	$interaction = substr($seq,0 ,$boxPos); 
    }

    return $interaction;
}


#######################
sub readFa
#reads given fa-file
#returns header and sequence
#######################
{
    my $file = shift;
    my ($seq,$head);
    my $fasta;

    open(FA, "<$file") or die $! ;

    while(<FA>){  
	$head =$_;
	if($head =~ /^>/){
	    $seq = uc(<FA>);
	    push(@{$fasta},$head.$seq);
	}
	else{
	    print $file," DOESN'T LOOK LIKE A FASTA FILE!!!!!\n";
	    exit(0);
	}
    }
    close(FA);
    
    return $fasta;
}

#############################################################################
sub writeDummyProfile
#write profile InteractionRegion_openen to profile folder
#needed by RNaplex
#############################################################################
{

    open(DUMMY,">".$opts{p}."/InteractionRegion_openen") or die print "WAS NOT ABLE TO CREATE DUMMY_PROFILE\nDO YOU HAVE PERMISSION TO WRITE TO FOLDER: ",$opts{p},"?\n";

    my @tmp = ("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA");
    print DUMMY "#opening energies\n";
    print DUMMY "#i\$\tl=1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\n";
    for(my $i =1;$i<21;$i++){

	$tmp[$i-1] = "0.000";
	print DUMMY $i,"\t",join("\t",@tmp),"\n";
    }

    close(DUMMY);
}
