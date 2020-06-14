#!/usr/bin/perl
use strict;
use List::Util qw(max);
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);


#Copyright (C) 2017-2019 University of Nagoya
#Contig_Hunter is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with Contig_Hunter.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################################################################################
#
#						Contig_Hunter version 0.2
#
#						Kazuma Uesaka
#						University of Nagoya
#						26 January 2017
#		
#						A Perl scripts to predict the pair of closest contig.
#
#	
#						Requirement:SMALT version 0.7.6 or higher (http://www.sanger.ac.uk/science/tools/smalt-0)
#
#		
#						Input:
#							 fastq file
#
#						Outnput:					
#							pair.txt to STDOUT
#						
#						Usage:
#						Contig_Hunter.pl (If you already made the path to this script).
#						perl Contig_Hunter.pl (If you dont made the path to this script).
#
#				The R1.fastq,R2.fastq and it's assembled contig.fasta should be included in the same folder,
#				
#
#
#
#	2017-01-26 version0.2
#	
#
##sleepをなんども入れているのは、ゲノムがシンプルでファイルが軽く、速く進みすぎるとなぜか出力がゼロになることが多発したため。
#version 1.4 maximum coverage support
#######################################################################################################################################


print "\n\n############################################################################################################################\n";
print "Program: Contig_Hunter\n";
print "version 1.0\n\n";
print "\nUsage:Contig_Hunter.pl <options>\n\n";
print "Input/output options:\n\n";
print "\t -1	R1.fastq (Required)\n";
print "\t -2	R2.fastq (Required)\n";
print "\t -f	contig/scaffold.fasta (Required)\n";
print "\t -o	Output file name (default results)\n";
print "\t -m	Run or skip mapping analysis (default 1)\n";
print "\t			1;runs\n";
print "\t			2;skip (mode for re-perform this analysis with different Matching Parameter.) \n";
print "\t -r	Average read size (bp) for calculating coverage (default 300)\n";

print "\nMatching Parameter options:\n\n";
print "\t -e	Sensitivity for call heterogous contig pair(default >= 2)\n";
print "\t -l	Minimum contig size (bp) to find pair (default >= 300)\n";
print "\t -s	Length (bp) to extract read from the end of contig (default <= 500)\n";
print "\t -c	Low threshould coverage to remove poorly mapped contig. If c >=1, mapping analysis perfromed again (default 0)\n";
print "\t -p	Maximum size to find heterogous mapped pair (default 1000)\n";


print "\nMapping Parameter options:\n\n";
print "\t -n	Number of threads for SMALT mapping (default 20)\n";
print "\t -h	Word length for SMALT mapping (default 14)\n";
print "\t -k	Skip step for SMALT mapping (default 8)\n";
print "\t -y	Identity threshold for SMALT mapping (default 150)\n";
print "\t -mp	Mismatch penalty for SMALT mapping (default -10)\n";
print "\t -go	Gap open penalty for SMALT mapping (default -10)\n";
print "\t -ge	Gap extension penalty for SMALT mapping (default -10)\n";
print "\t -force	Selct best hit pair and discard other pairs (default 0)\n";
print "\t			0;not run \n";
print "\t			1;run \n";

print "\nCircular contig options:\n\n";
print "\t -sl	Minumum size (bp) to call contig as circular (default 1000)\n";
print "\t -ml	Maximum size (bp) to call contig as circular (default 100000)\n";
print "\t -minc	Low threshould coverage to call contig as circular (default 0)\n";
print "\t -maxc	Maximum threshould coverage to call contig as circular (default 100000)\n";
print "############################################################################################################################\n\n";
my @now = localtime;print "\n# $now[2]:$now[1]:$now[0] ";
print "Starts Contig_Hunter\n";
print "Usage: Contig_Hunter\n";
system("sleep 3s");
my $fastq1 = "";
my $fastq2 = "";
my $fasta = "";
my $output = "results";
my $edge = "2";
my $minsize = "300";
my $Psize = "500";
my $word = "14";
my $skip = "8";
my $ithreshold = "150";
my $MismatchPenalty = "-10";
my $GapOpenPenalty = "-10";
my $GapExtentionPenalty = "-10";
my $CPU = "20";
my $check = "1";
my $readsize = 301;
my $lowcoverage = 0;
my $maxnumber = 1000;
my $mincirclesize = 1000;
my $maxcirclesize = 100000;
my $mincirclecov = 0;
my $maxcirclecov = 100000;
my $mode = 0;

#GetOptions に渡すオプション名の後ろに =s を付けると、オプションに文字列引数をとることができるようになる。=i を付けると整数引数をとれる（型指定はもちろんこのほかにも用意されている）。
GetOptions('1=s' => \$fastq1,'2=s' => \$fastq2, 'f=s' => \$fasta, 'o=s' => \$output,'m=f' => \$check,'r=i' => \$readsize, 'e=f' => \$edge, 'c=f' => \$lowcoverage, 'l=f' => \$minsize, 's=f' => \$Psize, 'h=i' => \$word, 'k=i' => \$skip, 'y=f' => \$ithreshold, 'mp=f' => \$MismatchPenalty, 'go=f' => \$GapOpenPenalty, 'ge=f' => \$GapExtentionPenalty, 'n=i' => \$CPU, 'p=i' => \$maxnumber, 'sl=i' => \$mincirclesize, 'ml=i' => \$maxcirclesize, 'minc=f' => \$mincirclecov, 'force=f' => \$mode, 'maxc=f' => \$maxcirclecov);

die "\ninput fastq file must be necesarry !\n\n\n" if($fastq1 eq "" or $fastq2 eq "");
die "\ninput 1.fastq file must be necesarry !\n\n\n" if($fastq1 eq "");
die "\ninput 2.fastq file must be necesarry !\n\n\n" if($fastq2 eq "");
die "\ninput fasta file must be necesarry !\n\n\n" if($fasta eq "");
open LOG, ">$output\_log" or die "cant open log file!\n";
print LOG <<HERE;#ヒアドキュメント機能でプリント
Input/output options:
-1	$fastq1
-2	$fastq2
-f	$fasta
-o	$output
-m	$check
-r	$readsize

Matching Parameter options:
-e	$edge
-l	$minsize
-s	$Psize
-c	$lowcoverage
-p	$maxnumber

Mapping Parameter options:
-n	$CPU
-h	$word
-k	$skip
-y	$ithreshold
-mp	$MismatchPenalty
-go	$GapOpenPenalty
-ge	$GapExtentionPenalty

Circular contig options:
-sl	$mincirclesize
-ml	$maxcirclesize
-minc	$mincirclecov
-maxc	$maxcirclecov
HERE

system("mkdir temporary");
system("mkdir $output\_R_ggplot2");
my $checker = 0;
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Processing fasta\n";
&convert1;
system("sleep 3s");
&convert2("tempfile");

my @now = localtime;print "# $now[2]:$now[1]:$now[0] Start mapping\n";
system("smalt index -k $word -s $skip hp temporary/contig_all.fasta");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Skip mapping\n" unless($check == 1);
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Make R1.sam\n" if($check == 1);
system("smalt map -n $CPU -d -1 -y $ithreshold -S subst=$MismatchPenalty,gapopen=$GapOpenPenalty,gapext=$GapExtentionPenalty -o temporary/R1.sam hp $fastq1") if($check == 1);
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Make R2.sam\n" if($check == 1);
system("smalt map -n $CPU -d -1 -y $ithreshold -S subst=$MismatchPenalty,gapopen=$GapOpenPenalty,gapext=$GapExtentionPenalty -o temporary/R2.sam hp $fastq2") if($check == 1);
system("mv *smi *sma temporary/");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Count coverage\n";
&coverage;
system("sleep 3s");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Extract R1 read end\n";
&extract;
system("sleep 3s");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Remove low coverage contig\n" if($lowcoverage > 0);
&remove if($lowcoverage > 0);
system("sleep 3s");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Extract heterougous mapped pair-reads\n";
&pair1;
system("sleep 3s");
&pair2;
system("sleep 3s");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Summarizing data...\n";
&pair3;
system("sleep 3s");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Summarizing data...\n";
&pair4;
system("sleep 3s");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Print the size of contig\n";
&length;
system("sleep 5s");
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Find pair of different contig end\n";
&node_edge;
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Remove low coverage contig information\n" if($lowcoverage > 0);
&remove2 if($lowcoverage > 0);
my @now = localtime;print "# $now[2]:$now[1]:$now[0] Select best hit edge\n" if($mode == 1);
&besthit if($mode == 1);
&graph;
system("R --vanilla --slave  --args results_R_ggplot2/coverage.txt results_R_ggplot2/length.txt results_R_ggplot2/GC.txt length_coverage_GC.pdf < plot.R ");#グラフ描画
system("mv temporary/Low_coverage.txt temporary/Low_coverage_contig_list.txt");#次回解析でエラーを出さないよう排除する。
print "\nRessults are saved as file name \"$output\_pair.txt\"\, \"$output\_length.txt\", \"$output\_length.txt\",and \"$output\_circular_pair.txt\", respectively.\n\n";
exit;
####################################################################################################################################
### SUBROUTINES ###

sub convert1 {
	open INPUT1, "<$fasta" or die "cant open txt perl1 file!\n";
	open INPUT2, "<$fasta" or die "cant open txt perl1 file!\n";
	open (OUT1, '>temporary/tempfile');
	my $file = "";
	my $line2 = <INPUT2>;
	my $count = 1;
	while (my $line = <INPUT1>) {
	my $line2 = <INPUT2>;#1行先読みファイルの入力
	chomp($line);#改行を除く
	print OUT1 "\>$count\t" if($line =~ "\>");#先頭行を出力
	$count++ if($line =~ "\>");#
	next if($line =~ "\>");
	print OUT1 "$line";
	print OUT1 "\n"  if($line2 =~ "\>");#1行先読みファイルを識別に利用している
	}
	print OUT1 "\n";
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub convert2 {
my $file = $_[0];#サブルーチンの引数は$_[0]に収納されている。
open INPUT3, "<temporary/$file" or die "cant find convert2 file!\n";
open (OUT2, '>temporary/contig_all.fasta');
open (OUT3, '>temporary/contig_distribution.vcf');

open (OUT4, '>temporary/.number');
open OUT5, ">$output\_R_ggplot2/GC.txt";
open OUT6, ">$output\_GC.txt";
my $Linenumber = 0;

print OUT5 "GC\n";
print OUT6 "node\tGC\n";
while (my $line = <INPUT3>) {
	chomp($line);
	my @array = split(/\t/, $line);
	my $size = length($array[1]);#まず長さを調べる。
	
	$array[0] =~ s/\>//;#vcfファイルには>は出力したくないので消しておく。
	next if ($size < $minsize);#default 300bp以上
	print OUT3 "$array[0]\t$size\n";#contig_distribution.vcfファイル出力、1カラム目は名前、2カラム目は長さ
	print OUT2 "\>$array[0]\n$array[1]\n";
	
	$Linenumber++;
	#次にGC含量計算
	my $Acount = 0;
	my $Tcount = 0;
	my $Gcount = 0;
	my $Ccount = 0;
	$Acount++ while($array[1] =~ m/A/g);
	$Tcount++ while($array[1] =~ m/T/g);
	$Gcount++ while($array[1] =~ m/G/g);
	$Ccount++ while($array[1] =~ m/C/g);
	my $GC = ($Gcount + $Ccount) / $size;
	my $GCratio = int $GC * 1000;#一回10倍して
	$GCratio = $GCratio / 10;#1/10にする。
	print OUT5 "$GCratio\n";
	
	print OUT6 "NODE\_$array[0]";
	print OUT6 "L\t$GCratio\n";
	print OUT6 "NODE\_$array[0]";
	print OUT6 "R\t$GCratio\n";
}#while終了
print OUT4 "$Linenumber";
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub coverage {
#contigの数だけループを繰り返して、全contigに対するカバレッジを順番に計算していく
open FA, "<temporary/contig_all.fasta" or die&&exit;
open NUM, "<temporary/.number" or die&&exit;
open OUT4, ">$output\_coverage.txt";
open OUT5, ">temporary/.cov";
open OUT5L, ">temporary/Low_coverage.txt";
open OUT6, ">$output\_R_ggplot2/coverage.txt";
my $size = <NUM>;
my @box = "";
for ( my $a = 1; $a <= $size; $a++ ){
	$box[$a] = 0;
}

open SAM1, "<temporary/R1.sam" if($checker == 0);
open SAM1, "<temporary/re-R1.sam" if($checker == 1);
while (my $seq1 = <SAM1>){#while開始
	next if($seq1 =~ /^\@/);
	my @array1 = split(/\t/, $seq1);
	$box[$array1[2]]++;
}

open SAM2, "<temporary/R2.sam" if($checker == 0);
open SAM2, "<temporary/re-R2.sam" if($checker == 1);
while (my $seq2 = <SAM2>){#while開始
	next if($seq2 =~ /^\@/);
	my @array1 = split(/\t/, $seq2);
	$box[$array1[2]]++;
}


print OUT4 "node\tcoverage\n";
print OUT6 "coverage\n";
while (my $name = <FA>){
	chomp($name);#1行目は名前
	$name =~ s/\>//;#>は消しておく。
	my $sequence = <FA>;chomp($name);#2行目は配列
	my $seqsize = length($sequence);
	my $avg = 0;
	$avg = int ($box[$name] / $seqsize * $readsize * 10);#10倍してからintで整数部分だけ返す
	$avg = $avg / 10;#10倍を解除して小数点1桁まで求める
	print OUT4 "NODE\_$name";
	print OUT4 "L\t$avg\n";
	print OUT4 "NODE\_$name";
	print OUT4 "R\t$avg\n";
	print OUT5 "$name\t$avg\n";
	print OUT5L "$name\t$avg\n" if($lowcoverage > 0 && $avg < $lowcoverage);#low coverage contigの出力
	print OUT6 "$avg\n";
}
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub extract {

#contig末端の指定サイズ以内にmapされているリードを抽出。
#step1のcontig_distribution.vcfを利用する
open INPUT1, "<temporary/contig_distribution.vcf" or die "cant open file!\n";
#coting_length.vcfを配列に収納

	my @name = (); my @length = ();my @right = ();my @left = ();my $i = 1;
	while (my $line = <INPUT1>) {
		chomp($line);
		my @array = split(/\t/, $line);
		$name[$i] = $array[0];
		$length[$i] = $array[1];
		$right[$i] = $array[1] - $Psize;#下流側の末端500bp以内
		$left[$i] = $Psize;#上流側の末端500bp以内
		$right[$i] = $array[1] if($array[1] < $Psize);#指定長（default 500bp）以下の短いcontigの場合は全長
		$left[$i] = $array[1] if($array[1] < $Psize);#上と同じく
		$i++;
	}

	#contigの末端のリードかどうかチェックする。
	open INPUT2, "<temporary/R1.sam";
	open INPUT3, "<temporary/R2.sam";
	open OUT5, ">temporary/R1out.txt" or die "cant open file!\n";

	while (my $line = <INPUT2>) {
		next if($line =~ /^\@/);
		chomp($line);
		my @array = split(/\t/, $line);
		next if($array[2] =~ /\*/);
		next unless($array[2] eq $name[$array[2]]);#文字列の完全一致は=~ではダメでeqを使う。=~は部分一致。
		print OUT5 "$array[0]\t$array[2]R\n" if($array[3] >= $right[$array[2]] && $array[1] == 0);
		print OUT5 "$array[0]\t$array[2]R\n" if($array[3] >= $right[$array[2]] && $array[1] == 256);#256は複数箇所
		print OUT5 "$array[0]\t$array[2]L\n" if($array[3] < $left[$array[2]] && $array[1] == 16);
		print OUT5 "$array[0]\t$array[2]L\n" if($array[3] < $left[$array[2]] && $array[1] == 272);
	}

my @now = localtime;print "# $now[2]:$now[1]:$now[0] Extract R2 read end\n";
	open OUT6, ">temporary/R2out.txt" or die "cant open file!\n";
	while (my $line = <INPUT3>) {
		next if($line =~ /^\@/);
		chomp($line);
		my @array = split(/\t/, $line);
		next if($array[2] =~ /\*/);
		next unless($array[2] eq $name[$array[2]]);#文字列の完全一致は=~ではダメでeqを使う。=~は部分一致。
		print OUT6 "$array[0]\t$array[2]R\n" if($array[3] >= $right[$array[2]] && $array[1] == 0);
		print OUT6 "$array[0]\t$array[2]R\n" if($array[3] >= $right[$array[2]] && $array[1] == 256);
		print OUT6 "$array[0]\t$array[2]L\n" if($array[3] < $left[$array[2]] && $array[1] == 16);
		print OUT6 "$array[0]\t$array[2]L\n" if($array[3] < $left[$array[2]] && $array[1] == 272);
	}
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub remove {
system("cp temporary/R1out.txt temporary/R1");
system("cp temporary/R2out.txt temporary/R2");
system("cp $output\_R_ggplot2/GC.txt $output\_R_ggplot2/.GC1");
system("cp $output\_GC.txt temporary/GC2");
system("sleep 5s");
open INPUT1, "<temporary/R1" or die "cant open R2 file!\n";
open INPUT2, "<temporary/R2" or die "cant open R2 file!\n";
open INPUT3, "<temporary/Low_coverage.txt";
open INPUT4, "<$output\_R_ggplot2/.GC1";
open INPUT5, "<temporary/GC2";
open OUT1, ">temporary/R1out.txt" or die "cant save R1out.txt file!\n";
open OUT2, ">temporary/R2out.txt" or die "cant save R2out.txt file!\n";
open OUT4, ">$output\_R_ggplot2/GC.txt" or die "cant save _R_ggplot2/GC.txt file!\n";
open OUT5, ">$output\_GC.txt" or die "cant save GC.txt file!\n";
open ZZZ1, ">temporary/R1remove.txt" or die "cant save R1remove.txt file!\n";
open ZZZ2, ">temporary/R1remove.txt" or die "cant save R1remove.txt file!\n";
my @name = ();my $i = 1;
while (my $line = <INPUT3>) {
	chomp($line);
	my @array = split(/\t/, $line);
	$name[$i] = $array[0];
	$i++;
}

while (my $line1 = <INPUT1>) {
	my @array1 = split(/\t/, $line1);
	$array1[1] =~ s/L//g;#Lを消す
	$array1[1] =~ s/R//g;#Rを消す
	my $match = 0;
	for(my $a = 1;$a <= $i;$a++){
		$match++ if($array1[1] == $name[$a]);
	}
	print ZZZ1 "$match\t$line1" if($match >= 1);
	print OUT1 "$line1" if($match == 0);#low coverage contigは出力されない
}
my @now = localtime;print "# $now[2]:$now[1]:$now[0] R1 low coverage contig removed\n";
while (my $line2 = <INPUT2>) {
	my @array2 = split(/\t/, $line2);
	$array2[1] =~ s/L//g;#Lを消す
	$array2[1] =~ s/R//g;#Rを消す
	my $match = 0;
	for(my $a = 1;$a <= $i;$a++){
		$match++ if($array2[1] == $name[$a]);
	}
	print ZZZ2 "$match\t$line2" if($match >= 1);
	print OUT2 "$line2" if($match == 0);#low coverage contigは出力されない
}
my @now = localtime;print "# $now[2]:$now[1]:$now[0] R2 low coverage contig removed\n";
my $line4 = <INPUT4>;my $c = 1;
print OUT4 "GC\n";
while (my $line4 = <INPUT4>) {
	my $match = 0;
	for(my $a = 1;$a <= $i;$a++){
		$match++ if($c == $name[$a]);
	}
	print OUT4 "$line4" if($match == 0);
	$c++;
}

my $line5L = <INPUT5>;
print OUT5 "node\tGC\n";
while (my $line5L = <INPUT5>) {
	my $line5R = <INPUT5>;
	$line5L =~ s/(NODE_)//g;
	$line5R =~ s/(NODE_)//g;
	my @array5L = split(/\t/, $line5L);
	my @array5R = split(/\t/, $line5R);
	$array5L[0] =~ s/L//g;#Lを消す
	$array5R[0] =~ s/R//g;#Rを消す
	my $match = 0;
	for(my $a = 1;$a <= $i;$a++){
		$match++ if($array5L[0] == $name[$a]);
	}
	print OUT5 "NODE\_$line5L" if($match == 0);
	print OUT5 "NODE\_$line5R" if($match == 0);
}

}#サブルーチン終了
#-----------------------------------------------------------------------------------------------------------------------------------
sub pair1 {
#R2outを1/20ずつに分解出力
open INPUT1, "<temporary/R2out.txt" or die "cant open pair1 R2out.txt1 file!\n";
my $i = 0;
while (my $line = <INPUT1>) {#行数を調べる
	$i++;
}
my $divide = $i / 20;

open INPUT2, "<temporary/R2out.txt" or die "cant open pair1 R2out.txt2 file!\n";
for(my $l = 1;$l <= 20;$l++){#20に分割する
	my $name = "R2_".$l;#name
	open (OUT, ">temporary/$name");#1/20ごとに別ファイル出力
	my $ii = 0;
	while (my $line2 = <INPUT2>) {
		print OUT "$line2";
		$ii++;
		last if($ii>$divide);
	}
}
system("cp /Users/user/perl_script/step4-2.pl ./temporary");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_1 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_2 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_3 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_4 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_5 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_6 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_7 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_8 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_9 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_10 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_11 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_12 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_13 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_14 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_15 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_16 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_17 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_18 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_19 &");
system("perl /Volumes/macbookpro_2015_backup/Users/kazumaxneo/perl_script/step4-2.pl R2_20");
system("sleep 5s");

system("cp temporary/NODE1 temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE2 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE3 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE4 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE5 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE6 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE7 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE8 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE9 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE10 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE11 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE12 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE13 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE14 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE15 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE16 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE17 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE18 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE19 >> temporary/another_contig_paired2.txt");system("sleep 2s");
system("cat temporary/NODE20 >> temporary/another_contig_paired2.txt");system("sleep 5s");
system("rm temporary/R2_* temporary/NODE*");system("sleep 2s");
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub pair2 {
#今ままでprocess5ではnode1Rxnode23Rなどを、node1Rxnode23R 10、node23Rxnode1R 13、などと二回別のものとしてカウントしていた。それを修正する。
#process5でまとまったことを確認。
#2 後半にunpairR2を出力するスクリプトを追加

open INPUT2, "<temporary/another_contig_paired2.txt" or die "cant open sam file!\n";
open OUT2, ">temporary/temp1";
while (my $line = <INPUT2>) {
	chomp($line);
	my @array = split(/\t/, $line);
	my $charaL = substr($array[0],-1);
	my $charaR = substr($array[1],-1);
	
	if($charaL =~ /R/ && $charaR =~ /L/){
		print OUT2 "$array[1]\t$array[0]\n";
	}else{
		print OUT2 "$array[0]\t$array[1]\n";
	}
}
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub pair3 {
#順番を入れ替えるだけ
open INPUT3, "<temporary/temp1" or die "cant open sam file!\n";
open OUT3, ">temporary/out.txt";
while (my $line3 = <INPUT3>) {
	chomp($line3);
	my @array = split(/\t/, $line3);
	my $left = $array[0];
	my $right = $array[1];
	$left =~ s/L//g;#Lを消す
	$right =~ s/L//g;#Rを消す
	if($left > $right){
		print OUT3 "NODE_$array[1]\tNODE_$array[0]\n";
	}else{
		print OUT3 "NODE_$array[0]\tNODE_$array[1]\n";
	}
}
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub pair4 {
#ネットで助けてもらって改善
open INPUT, "<temporary/out.txt" or die "cant open sam file!\n";

my @boxLR = ();
my @box = ();
while (my $l = <INPUT>) {
	chomp($l);
	$l =~ s/(NODE)//g;
	my @div = split(/\_/, $l);
	my @div2 = split(/\t/, $l);
	#contig9_L	contig7_R
	unless($div2[0] =~ /\_/){#nameにLとRの区別がないものをここで配列収納。
		push(@box,$l);
	}else{	#nameにLとRの区別があるものを次に配列収納。
		push(@boxLR,$l) if($div[1] =~ "L");
		push(@boxLR,$l) if($div[1] =~ "R");
	}
}
unshift @boxLR, @box;#前の配列に連結

open IA, ">temporary/summarized.txt" or die "cant open sam file!\n";
my $ref;
foreach (@boxLR){
chomp;
my @t = split /\t/;
$ref->{$t[0]}->{$t[1]}++ ;#$refの最後に++が実行され、%refのバリューに入る。++はその組み合わせが複数出てきた時に実行されるカウンタとなっている。例えばcontig10_L x contig10_Rの組み合わせが見つかるたびに1数字がカウントされる。
}
for my $v1 (sort keys %$ref) {#ハッシュのハッシュである%refのキー（要素）は$t[0]なので、$t[0]の若い順でソートされる。
for my $v2 (sort keys %{$ref->{$v1}}) {#%refの中のハッシュである{$t[0]}->{$t[1]}のキー（要素）は$t[1]なので、さらに$t[1]の若い順でソートされる。
print IA "NODE$v1\tNODE$v2\t$ref->{$v1}->{$v2}\n";
}
}

}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub length {
open INPUT1, "<temporary/contig_distribution.vcf" or die;
open INPUT2, "<temporary/Low_coverage.txt";
open OUTPUT1, ">$output\_length.txt" or die;
open OUTPUT2, ">$output\_R_ggplot2/length.txt" or die;
open OUTPUT3, ">$output\_R_ggplot2/name.txt" or die;
my $i = 1;
print OUTPUT1 "node\tlength\n";
print OUTPUT2 "length\n";
print OUTPUT3 "name\n";
#1変数準備用whileループ
my @name = ();my $i = 1;
while (my $line = <INPUT2>) {
	chomp($line);
	my @array = split(/\t/, $line);
	$name[$i] = $array[0];
	$i++;
}

while (my $in = <INPUT1>){
	chomp($in);
	my @array = split(/\t/,$in);
	my $match = 0;
	for(my $a = 1;$a <= $i;$a++){
		$match++ if($array[0] == $name[$a]);
	}
	next unless($match == 0);#low coverage filterがある場合、low coverage相当のcontig情報は出力しない。
	print OUTPUT1 "NODE\_$array[0]L\t$array[1]\n";
	print OUTPUT1 "NODE\_$array[0]R\t$array[1]\n";
	print OUTPUT2 "$array[1]\n";
	print OUTPUT3 "NODE\_$array[0]\n";
	}#2while終了
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub node_edge {
open INPUT, "<temporary/summarized.txt" or die;
open INPUT2, "<temporary/.cov" or die;
open INPUT3, "<temporary/contig_distribution.vcf" or die "cant open file!\n";
open INPUT4, "<temporary/Low_coverage.txt";
open NUM, "<temporary/.number" or die&&exit;
open OUTPUT2, ">$output\_pair.txt" or die;
open OUTPUT3, ">$output\_circular_pair.txt" or die;

print OUTPUT2 "node1\tedge\tnode2\n";
print OUTPUT3 "node\tedge\tlength\tcoverage\n";
my $contigsize = <NUM>;
my @name = (); my @length = ();my $i = 0;
while (my $line = <INPUT3>) {
	chomp($line);
	my @array = split(/\t/, $line);
	$name[$i] = $array[0];
	$length[$i] = $array[1];
	$i++;
}

my @name2 = (); my @cov = ();my $ii = 0;
while (my $line = <INPUT2>) {
	chomp($line);
	my @array = split(/\t/, $line);
	$name2[$ii] = $array[0];
	$cov[$ii] = $array[1];
	$ii++;
}

my @name4 = ();my $iii = 1;
while (my $line4 = <INPUT4>) {
	chomp($line4);
	my @array4 = split(/\t/, $line4);
	$name4[$iii] = $array4[0];
	$iii++;
}


#1変数準備用whileループ
while (my $in = <INPUT>){
	chomp($in);
	$in =~ s/(NODE_)//g;
	my @array = split(/\t/,$in);
	next unless ($array[2] >= $edge);#ペアリードの個数。default2以上。
	print OUTPUT2 "NODE\_$array[0]\t$array[2]\tNODE\_$array[1]\n";
	$array[0] =~ s/L//g;#Lを消す
	$array[0] =~ s/R//g;#Rを消す
	$array[1] =~ s/L//g;#Lを消す
	$array[1] =~ s/R//g;#Rを消す
	next unless($array[0] eq $array[1]);#環状になるcontigかどうか判定
	for(my $a = 0;$a <= $i;$a++){
		next unless($array[0] eq $name[$a]);
		for(my $b = 0;$b <= $ii;$b++){
			print OUTPUT3 "NODE\_$array[0]\t$array[2]\t$length[$a]\t$cov[$b]\n" if($array[0] eq $name[$b] && $length[$a] >= $mincirclesize && $length[$a] <= $maxcirclesize && $cov[$b] >= $mincirclecov && $cov[$b] <= $maxcirclecov);
		}
	}
}#2while終了


#自分自身(LとR)との繋がりのデータを出力。(e.g., NODE_1L	10000	NODE_1R)
for(my $c = 1;$c <= $contigsize;$c++){
	my $match = 0;
	for(my $d = 1;$d <= $iii;$d++){
		$match++ if($c eq $name4[$d]);
	}
	print OUTPUT2 "NODE\_$c" if($match == 0);
	print OUTPUT2 "L\t10000\tNODE\_$c" if($match == 0);
	print OUTPUT2 "\R\n" if($match == 0);
}
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub remove2 {
system("cp $output\_R_ggplot2/coverage.txt $output\_R_ggplot2/.cov1");
system("cp $output\_coverage.txt temporary/.cov2");
system("sleep 3s");
open INPUT1, "<temporary/Low_coverage.txt";
open INPUT2, "<$output\_R_ggplot2/.cov1";
open INPUT3, "<temporary/.cov2";
open OUT1, ">$output\_R_ggplot2/coverage.txt" or die "cant save R_ggplot2/coverage.txt file!\n";
open OUT2, ">$output\_coverage.txt" or die "cant save coverage.txt file!\n";

my @name = ();my $i = 1;
while (my $line = <INPUT1>) {
	chomp($line);
	my @array = split(/\t/, $line);
	$name[$i] = $array[0];
	$i++;
}

my $line2 = <INPUT2>;my $b = 1;
print OUT1 "coverage\n";
while (my $line2 = <INPUT2>) {
	my $match = 0;
	for(my $a = 1;$a <= $i;$a++){
		$match++ if($b == $name[$a]);
	}
	print OUT1 "$line2" if($match == 0);
	$b++;
}

my $line3L = <INPUT3>;
print OUT2 "node\tcoverage\n";
while (my $line3L = <INPUT3>) {
	my $line3R = <INPUT3>;
	my @array3L = split(/\t/, $line3L);
	my @array3R = split(/\t/, $line3R);
	$array3L[0] =~ s/(NODE_)//g;
	$array3R[0] =~ s/(NODE_)//g;
	$array3L[0] =~ s/L//g;#Lを消す
	$array3R[0] =~ s/R//g;#Rを消す
	my $match = 0;
	for(my $a = 1;$a <= $i;$a++){
		$match++ if($array3L[0] == $name[$a]);
	}
	print OUT2 "$line3L" if($match == 0);
	print OUT2 "$line3R" if($match == 0);
}

}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub besthit {
#results_pair.txtを使う。ベストヒットだけに絞る。
system("sleep 5s");
open INPUT1, "<$output\_pair.txt" or die "cant open file!\n";
open INPUT2, "<$output\_pair.txt" or die "cant open file!\n";
open INDEX, ">index" or die "cant open file!\n";

#重複度をインデックス化
my $i = 1;my $line1 = <INPUT1>;my $line1 = <INPUT1>;
my $line2 = <INPUT2>;
while (my $line2 = <INPUT2>) {
		my $line1 = <INPUT1>;
		chomp($line1);chomp($line2);
		my @array1 = split(/\t/, $line1);
		my @array2 = split(/\t/, $line2);
		print INDEX "$i\n" unless($array1[0] eq $array2[0]);
		$i = 1 unless($array1[0] eq $array2[0]);
		$i++ if($array1[0] eq $array2[0]);
}
system("sleep 3s");

open INDEX2, "<index" or die "cant open file!\n";
open INPUT3, "<$output\_pair.txt" or die "cant open file!\n";
open OUT, ">$output\_best_hit_pair.txt" or die "cant open file!\n";
my $line3 = <INPUT3>;

print OUT "node1\tedge\tnode2\n";
while (my $line3 = <INPUT3>) {
	my $index = <INDEX2>;
	chomp($line3);chomp($index);
	my @array3 = split(/\t/, $line3);
	for(my $a = 1;$a <= $index;$a++){
		print OUT "$line3\n" if($index == 1);
		last if($index == 1);
		my @array4 = "";
		if ($a == $index){
			$array4[1] = 0;
		}else{
			my $next = <INPUT3>;
			chomp($next);
			@array4 = split(/\t/, $next);
		}
		if($array4[1] > $array3[1]){
			$array3[1] = $array4[1];
			$array3[2] = $array4[2];
		}
	}
	print OUT "$array3[0]\t$array3[1]\t$array3[2]\n" unless($index == 1);
}
}#サブルーチン終了


#-----------------------------------------------------------------------------------------------------------------------------------
sub graph {
open OUT, "> plot.R" or die "cant save file!\n";
my $text = <<"EOS";
library(ggplot2)
library(scales)
args <- commandArgs(trailingOnly = T)
( x <- read.table(args[1], header=T) )
( y <- read.table(args[2], header=T) ) #x軸も同様に入力
( z <- read.table(args[3], header=T) ) #3番目の要素
df <- data.frame(x=x, y=y,z)
g <- ggplot(df,aes(x = x,y = y,z))
g <- g + xlab("coverage")    # x 軸ラベル
g <- g + ylab("length")    # y 軸ラベル
g + geom_point(aes(colour=GC,size = coverage),alpha = 0.5) + scale_size_area(name = "coverage", max_size = 8) + scale_x_log10(labels=trans_format("log10",math_format(10^.x))) + scale_y_log10(labels=trans_format("log10",math_format(10^.x))) + scale_colour_gradientn(colours=rainbow(3)) + theme(axis.title.x = element_text(size=15)) + theme(axis.title.y = element_text(size=15)) + theme(axis.text.x = element_text(size=15)) + theme(axis.text.y = element_text(size=15)) + coord_cartesian(xlim=c(0.1,5000),ylim=c(100,2000000))
ggsave(file=args[4])#セーブ
dev.off()
EOS
print OUT "$text\n";
}
####################################################################################################################################