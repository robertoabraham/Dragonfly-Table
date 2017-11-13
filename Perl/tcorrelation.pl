#!/opt/local/bin/perl

use Getopt::Long;
use Pod::Usage;

# Parse command-line options
my $help = 0;
my $man = 0;

$help = 1 if $#ARGV == -1;
$result = GetOptions(
    "help|?" => \$help, 
     man=> \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus =>0, -verbose =>2) if $man;

# Load the columns
my $keyX = $ARGV[0];
my $keyY = $ARGV[1];

# Load the data table
my $xmax = -1e100;
my $xmin =  1e100;
my $ymax = -1e100;
my $ymin =  1e100;
my $colnum=0;
my $count=0;

while(<STDIN>)
{
    my $line = $_;
    if ($line =~ /^#/) {

        next if $line =~ /^$\!/;

        @fields = split;
        $keyword = $fields[2]; 
        
        # Define column numbers to correspond to SExtractor keywords
        # so for example: $col{"FLUX_ISO"}=3 etc.
        $col{$keyword}=$colnum++;

        next;
    }
    chomp($line);
    $line =~ s/^\s+//;
    @column = split(' ',$line);

    $x[$count] = $column[$col{$keyX}];
    $x[$count] = log($x[$count])/log(10) if ($xlog);
    $xmin = $x[$count] if ($x[$count] !~ /n/ && $x[$count] < $xmin);
    $xmax = $x[$count] if ($x[$count] > $xmax);

    $y[$count] = $column[$col{$keyY}];
    $y[$count] -= $skylevel if $skylevel;
    if ($ylog){
        next if ($y[$count] <= 0);
        $y[$count] = log($y[$count])/log(10);
    }
    $ymin = $y[$count] if ($y[$count] !~ /n/ && $y[$count] < $ymin);
    $ymax = $y[$count] if ($y[$count] > $ymax);

    $count++;
}


# Mean values
$mx = 0;
$my = 0;
for ($i=0;$i<$count;$i++) {
    $mx += $x[$i];
    $my += $y[$i];
}
$mx /= $count;
$my /= $count;

# Correlation coefficients (un-normalized second moments)
$sxx = 0;
$sxy = 0;
$syy = 0;
for ($i=0;$i<$count;$i++) {
    $sxx += ($x[$i]-$mx)*($x[$i]-$mx);
    $syy += ($y[$i]-$my)*($y[$i]-$my);
    $sxy += ($x[$i]-$mx)*($y[$i]-$my);
}

# Variances and covariances (normalized correlation coefficients)
$vx  = $sxx/$count;
$vy  = $syy/$count;
$cov = $sxy/$count;

# The covariance matrix is:
#  [vx  cov]
#  [cov  vy] 
#
#  If T is the trace and D is the determinant, the eigenvalues are:
#
#  L1 = T/2 + ((T**2)/4-D)**(1/2)
#  L2 = T/2 - ((T**2)/4-D)**(1/2)

$T = ($vx + $vy);
$D = ($vx * $vy - $cov * cov);
$L1 = $T/2.0 + (($T**2)/4.0-$D)**0.5;
$L2 = $T/2.0 - (($T**2)/4.0-$D)**0.5;

printf("- Moments -\n");
printf("Mean (column 1):     %f\n",$mx); 
printf("Variance (column 1): %f\n",$vx); 
printf("\n");
printf("Mean (column 2):     %f\n",$my); 
printf("Variance (column 2): %f\n",$vy); 
printf("\n");
printf("Correlation:         %f\n",$count * $cov);
printf("Covariance:          %f\n",$cov);
printf("\n");
printf("- Covariance matrix -\n");
printf("[ %f     %f ]\n",$vx,$cov);
printf("[ %f     %f ]\n",$cov,$vy);
printf("\n");
printf("Trace:           %f\n",$T);
printf("Determinant:     %f\n",$D);
printf("Eigenvalue1:     %f\n",$L1);
printf("Eigenvalue2:     %f\n",$L2);
printf("\n");
printf("- Ellipse parameters -\n");
printf("Axis1:           %f\n",sqrt($L1));
printf("Axis2:           %f\n",sqrt($L2));

#printf("2.35*sqrt(L1) = %5.2f    2.35*sqrt(L2) = %5.2f\n",2.35 * sqrt($L1),2.35 * sqrt($L2));
#printf("in arcsec     = %5.2f    in arcsec     = %5.2f\n",18.5*2.35 * sqrt($L1),18.5*2.35 * sqrt($L2));




__END__

=head1 NAME

tcorrelation - Compute correlation/covariance statistics for two columns of a data table

=head1 SYNOPSIS

tcorrelation [options] < data.txt 

=over 8

=item B<xaxis>

Column name for the X-axis.

=item B<yaxis>

Column name for the Y-axis.

=back

=head1 OPTIONS

=over 8

=item B<--help>

Prints a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<tcorrelation> computes the following:

    * mean
    * standard deviation
    * covariance matrix
    * eigenvalues of the covariance matrix
    * eigenvectors of the covariance matrix
    * parameters of the 1 sigma error ellipse

=cut
