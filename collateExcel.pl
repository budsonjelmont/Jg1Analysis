#!/usr/bin/perl -w

use diagnostics;
#use strict;
use Win32::OLE qw(in with);
use Win32::OLE::Const 'Microsoft Office .* Object Library';
use File::Find;

my $filesPath = $ARGV[0];
my $outfileName = $ARGV[1];
my $renameColumns = $ARGV[2];

#If $renameColumns = TRUE, then replace the column names in each file with the names below.
#Must be in the same order they appear in the original excel spreadsheets 
my @colNames = ("Protein name", "Phosphosite", "Database direction",
	"Assigned sequence", "Validation score", "Ascore", "MOWSE score", "Delta mass (ppm)",
	"Isolated mass", "Scan number", "Charge state", "UNIPROT accession",
	"HPRD accession", "All database identifiers containing peptide sequence", "Peak area") ;

#get list of Excel files to collate
my @files;
my $start_dir = $filesPath;  # top level dir to search
find( 
    sub { push @files, $File::Find::name unless -d; }, 
    $start_dir
);

#create outfile in same directory
my $excelOut = $filesPath."\\".$outfileName.".xls";

#create Excel object
my $Excel = Win32::OLE->GetActiveObject('Excel.Application')
        || Win32::OLE->new('Excel.Application', 'Quit');
#stop if errors
$Win32::OLE::Warn = 3;   
#turn off warnings/alerts
$Excel->{DisplayAlerts}=0;

#create new excel workbook for collated files...
my $outBook = $Excel->Workbooks->Add();
   $outBook->SaveAs($excelOut);
#...or open it if it's already been created
$outBook = $Excel->Workbooks->Open($excelOut); 

for(my $i=0; $i<=$#files; $i++){
	print "Now working on $files[$i]\n";
	# open Excel file
	my $inBook = $Excel->Workbooks->Open($files[$i]);

	# select worksheet, rename it, and fix column names
	my $inSheet = $inBook->Worksheets(1);
#	$files[$i] =~ /[\/]{1}([a-zA-Z0-9_\s]+)\.xls/;	#matches whole file name, but if filename is > ~30 chars, excel won't accept it
	$files[$i] =~ /Jurkat_PLCr_((WT|KO)_IAP_[0-9]+min_R[0-9]+)/;
	$inSheet->{Name} = $1;
	
	if($renameColumns eq TRUE){
		for(my $col = 0; $col <= $#colNames; $col++){
			# $inBook = $Excel->Workbooks->Open($files[$i]);
			# $inSheet = $inBook->Worksheets(1);
			$inSheet->Cells(1,$col+1)->{'Value'} = $colNames[$col];
		}
	}
	
	#copy the data into the collated workbook
	#NOTE: collated workbook needs to be reopened every time, not sure why. If it isn't reopened/object isn't redefined,
	#the code below fails, saying that $outBook isn't defined.
	$outBook = $Excel->Workbooks->Open($excelOut);
	my $outSheet = $outBook->Sheets( $outBook->Sheets->Count );
	$inSheet->Move( { After => $outSheet });
	
	$outBook-> Save();
	
	#Close workbooks
	$inBook = $Excel->Workbooks->Close();
}

#Delete the first sheet in the collated workbook, since it's blank
#$outBook = $Excel->Workbooks->Open($excelOut);
$outBook -> Worksheets("Sheet1") -> Delete();
#Save and close the collated workbook
$outBook = $Excel->Workbooks->Close();