    use Sub::Recursive;

    my $fac = recursive {
        my ($n) = @_;
        return 1 if $n < 1;
        return $n * $REC->($n - 1);
    };
    
    print $fac->(5)."\n";