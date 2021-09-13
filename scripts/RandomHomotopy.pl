use Benchmark qw(:all);

#Try to reduce a complex to a point trough collapses and "anticollapses" applied trough the application R_H.cc

print "This defines a function that tries to reduce a complex to a point trough collapses and anticollapses applied trough the application R_H.cc\n",
      "use like, e.g.,\n",
      '  Rand_Hom(n,max,SC,sphere);',"\n",
      "where n is the number of rounds, max the maximum number of anticollapses, SC the simplicial complex\n",
      "and sphere=1 if the complex is a sphere, 0 otherwise\n";
      

sub Rand_Hom($$$$) {



use application 'topaz';


      
      
my ($rounds,$max,$inS,$sphere)=@_;


my $media=0;
my $fail=0;

my $minima=$max;
my $maximum=0;


 for (my $j=0;$j<$rounds;$j++) {
 
my $q=$inS; 
my $p=collapses($q);

if ($sphere>0) {
#delete a random facet
my $num_fac=$p->N_FACETS;

my @Fac=$p->FACETS;

$q=deletion($p,$Fac[0][int(rand($num_fac))]);

$p=collapses($q);
}

my $n=$p->DIM;
my $i=0;
while (($n >1)and($i<$max)) {

$q=R_H($p, $n, 1);

if ($q->DIM>1) {
$p=collapses($q);

$n=$p->DIM;

$i++;
}

if ($minima>$i) {
$minima=$i;
}

if ($maximum<$i) {
$maximum=$i;
}

if ($i==$max) {
$fail+=1;
}
else {
$media+=$i;
}

#return $p;
}
if ($fail<$rounds) {
$media=$media/($rounds -$fail);
}

 print "Rounds:".$rounds."\n";
 print "Media:".$media."\n";
 print "Failures:".$fail."\n";
 print "Minimum:".$minima."\n";
 print "Maximum:".$maximum."\n";
}
