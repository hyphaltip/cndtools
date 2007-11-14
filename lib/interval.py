class Interval:
    def __init__(self, chrom, start, end, strand='+'):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

    def __lt__(self, other):
        return self.chrom < other.chrom \
               or (self.chrom == other.chrom and self.start < other.start)

    def __sub__(self, other):
        return Interval(self.chrom, other.end, self.start, self.strand)

    def __str__(self):
        return "%s:%d-%d%s" % (self.chrom, self.start, self.end, self.strand)

    def contains(self, other):
        return (self.chrom == other.chrom and
                self.start <= other.start and
                self.end >= other.end)
