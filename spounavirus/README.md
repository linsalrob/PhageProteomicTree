###Phage Proteomic Tree for the Spounavirus analyses

The file spounavirus.tsv is the spreadsheet with all the spounavirus IDs in. We prefered the RefSeq ID if possible, but otherwise we used the INSDC ID. To make a file with just the IDs we used this one line perl code:


```perl -ne '@a=split /\t/; map {$a[$_] =~ s/^\s+//; $a[$_] =~ s/\s+$//} (1,2); if ($a[2] eq "-") {print "$a[1]\n"} else {print "$a[2]\n"}' spounavirus.tsv > spounavirus_ids.txt```
