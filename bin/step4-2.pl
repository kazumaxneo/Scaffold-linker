#!/usr/bin/perl
use strict;

#step3の計算ファイルを利用する
#異種contigにマップされたリードペアの情報を出力（マップされたcontig名とその座標）
my $in = $ARGV[0];
open INPUT4, "<temporary/$in" or die "cant open file!\n";
#配列に収納
my @name = (); my @map = ();my @position = ();my $b = 0;
while (my $line = <INPUT4>) {
	chomp($line);
	my @array = split(/\t/, $line);
	$name[$b] = $array[0];
	$map[$b] = $array[1];
	$b++;
}
print "$b\n";
$in =~ s/(R2_)//g;#出力ファイル名に使う。
#異種contigのペアかどうかチェックする。
open INPUT5, "<temporary/R1out.txt" or die "cant open file!\n";
open OUT1, ">temporary/NODE$in" or die "cant open file!\n";

my $i = 0;
while (my $line = <INPUT5>) {
		
		chomp($line);
		my @array = split(/\t/, $line);
		for(my $d = 0;$d <= $b;$d++){
			if($array[0] eq $name[$d]){
				#printの順番はR1のマッピングされたcontig、R1のマッピングされたcontigの座標、R2のマッピングされたcontig、R2のマッピングされたcontigの座標
				print OUT1 "$array[1]\t$map[$d]\n" unless($array[2] eq $map[$d]);
			}
		}
		$i++;
		print "$i\n";
}


exit;


