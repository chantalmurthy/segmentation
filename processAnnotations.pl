use List::Util qw(min max);
use List::MoreUtils qw( minmax );
use Scalar::Util qw(looks_like_number);
use POSIX;

##############
sub point_in_polygon { #order of specification of polygon vertices matters.
    my ( $x, $y, @xy ) = @_;
    my $n = @xy / 2;                      # Number of points in polygon.
    my @i = map { 2 * $_ } 0 .. (@xy/2);  # The even indices of @xy.
    my @x = map { $xy[ $_ ]     } @i;     # Even indices: x-coordinates.
    my @y = map { $xy[ $_ + 1 ] } @i;     # Odd indices: y-coordinates.

    my ( $i, $j );                        # Indices.

    my $side = 0;                         # 0 = outside, 1 = inside.

    for ( $i = 0, $j = $n - 1 ; $i < $n; $j = $i++ ) {
        if (
            (

             # If the y is between the (y-) borders ...
             ( ( $y[ $i ] <= $y ) && ( $y < $y[ $j ] ) ) ||
             ( ( $y[ $j ] <= $y ) && ( $y < $y[ $i ] ) )
            )
            and
            # ...the (x,y) to infinity line crosses the edge
            # from the ith point to the jth point...
            ($x
             <
             ( $x[ $j ] - $x[ $i ] ) *
             ( $y - $y[ $i ] ) / ( $y[ $j ] - $y[ $i ] ) + $x[ $i ] )) {
          $side = not $side; # Jump the fence.
      }
    }
    return $side ? 1 : 0;
}

##############
#Output array with each entry storing a string corresponding to tile's membership to an ROI
$scale = 45;
$width =  $ARGV[0]; ##CHANGE TO COMMAND LINE INPUT
$height =  $ARGV[1];
$nTilesWidth = ceil($width /$scale);
$nTilesHeight = ceil($height/$scale);
@outputArray = ("") x $nTilesWidth*$nTilesHeight; 

$annotationName = $ARGV[2];
chomp($annotationName);
$FileInput = "C:/Users/Nameeta/Chantal/data/W12-1-1-D.2/annotations/".$annotationName.".xml";

$TemplateResult = open(ANNOTATION, $FileInput);
$numlines = 0;
@curPoly;
#while(($tline = <ANNOTATION>) && ($numlines < 208)) {
while($tline = <ANNOTATION>) {
	@ttokens = split("<Annotation ", $tline);
	@ttokens2 = split("</Vertices>", $tline);
	@ttokens3 = split("<Vertex ", $tline);	
	if(@ttokens > 1){
		my $st = $ttokens[1];
		my @ar = split("Name=", $st);
		my $meh = $ar[1];
		my @bleh = split(">", $meh);
		$name = $bleh[0];
		$curROI = $name;
		chomp($curROI);
	} elsif (@ttokens2 > 1) {
		my $arrSize = @curPoly;
		#print "arrSize: $arrSize\n";
		#print "@curPoly\n";
		my $n = @curPoly / 2;                      # Number of points in polygon.
		my @i = map { 2 * $_ } 0 .. (@curPoly/2);  # The even indices of @curPoly.
		my @x = map { $curPoly[ $_ ]     } @i;     # Even indices: x-coordinates.
		my @y = map { $curPoly[ $_ + 1 ] } @i;     # Odd indices: y-coordinates.		
		pop(@x);
		pop(@y);	#some junk added to end of these twins.
		my ($xmin, $xmax) = minmax @x;
		my ($ymin, $ymax) = minmax @y;			
		#print "xmin: $xmin, ymin: $ymin, xmax: $xmax, ymax: $ymax\n";
		#tiles are one indexed
		$val1 = ceil($xmin/$scale)+1;
		$val2 = floor($xmax/$scale)-1;
		$val3 = ceil($ymin/$scale)+1;
		$val4 = floor($ymax/$scale)-1;	
		#print "$val1, $val2, $val3, $val4\n";
		for ($x_tile_ind=(floor($xmin/$scale)); $x_tile_ind<(ceil($xmax/$scale)); $x_tile_ind++) {
			for ($y_tile_ind=(floor($ymin/$scale)); $y_tile_ind<(ceil($ymax/$scale)); $y_tile_ind++) {				
				my $a = point_in_polygon(($x_tile_ind)*$scale, ($y_tile_ind)*$scale, @curPoly);
				my $b = point_in_polygon(($x_tile_ind)*$scale, ($y_tile_ind+1)*$scale, @curPoly);
				my $c = point_in_polygon(($x_tile_ind+1)*$scale, ($y_tile_ind)*$scale, @curPoly);
				my $d = point_in_polygon(($x_tile_ind+1)*$scale, ($y_tile_ind+1)*$scale, @curPoly);
				#check if 4 points lie within polygon	
				#if yes, store current ROI type in x_tile_ind, y_tile_ind in output array
				#if ($a || $b || $c || $d) {
				#	$outputArray[($y_tile_ind-1)*$nTilesWidth + $x_tile_ind] = $curROI;
				#}			
				if (($a && $b) || ($b && $c) || ($c && $d) || ($a && $c) || ($a && $d) || ($b && $d)) {
					$outputArray[($y_tile_ind-1)*$nTilesWidth + $x_tile_ind] = $curROI;
				} ##CHANGE TO COMMAND LINE INPUT to control number of corners to be shared
				print "abcd: $a $b $c $d\n";
			}
		}
		@curPoly = ();			
	} elsif (@ttokens3 > 1) {
		my $st = @ttokens3[1];
		my @ar = split("X=", $st);
		my @ar2 = split("Y=", $ar[1]);
		my $x = $ar2[0];
		$x = eval $x;
		my @bleh = split("/>" ,$ar2[1]);
		my $y = $bleh[0];	
		$y = eval $y;
		push(@curPoly, $x);
		push(@curPoly, $y);	
	}
	$numlines++;
}
close(ANNOTATION);

$FileOutput = "C:/Users/Nameeta/Chantal/results/Annotations/".$annotationName.".csv";

open(OUTPUT, ">$FileOutput");
$tileind = 1;
while($tileind <= $nTilesWidth*$nTilesHeight) {
	$st = $outputArray[$tileind];
	if ($st ne "") {
		chomp($st);
		$st = substr($st, 1, length($st)-2);
	} else {
		$st = 0;
	}	
	my $y = ceil($tileind/$nTilesWidth);
	my $x = $tileind - $nTilesWidth*(floor($tileind/$nTilesWidth));
	if ($x eq 0) {
		$x = $nTilesWidth;
	}
	
	$mode = $ARGV[3];
	chomp($mode);
	if ($mode eq "2") { 
		$x = $x*$scale;
		$y = $y*$scale;
	}
	if ($st ne 0) {
		print OUTPUT "$x,$y,$st,0,0,0\n";
	}
	$tileind = $tileind + 1;
}
close(OUTPUT);
