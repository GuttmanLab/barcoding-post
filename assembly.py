from collections import OrderedDict

class Assembly(object):
    _chromsizes = None

    def init_offsets(self, resolution):
        count = 0
        self._offsets = OrderedDict()
        for (chrom, size) in self._chromsizes.iteritems():
            self._offsets[chrom] = count
            count += -(-size // resolution)

    def get_size(self, chrom):
        return self._chromsizes.get(chrom)

    def get_offset(self, chrom):
        return self._offsets.get(chrom)

class Mm9(Assembly):

    def __init__(self):
        self._chromsizes = OrderedDict([
            ('chr1', 197195432),
            ('chr2', 181748087),
            ('chr3', 159599783),
            ('chr4', 155630120),
            ('chr5', 152537259),
            ('chr6', 149517037),
            ('chr7', 152524553),
            ('chr8', 131738871),
            ('chr9', 124076172),
            ('chr10', 129993255),
            ('chr11', 121843856),
            ('chr12', 121257530),
            ('chr13', 120284312),
            ('chr14', 125194864),
            ('chr15', 103494974),
            ('chr16', 98319150),
            ('chr17', 95272651),
            ('chr18', 90772031),
            ('chr19', 61342430),
            ('chrX', 166650296)])

    @classmethod
    def is_named(cls, name):
        return name == "mm9"

class Hg19(Assembly):

    def __init__(self):
        self._chromsizes = OrderedDict([
            ('chr1', 249250621),
            ('chr2', 243199373),
            ('chr3', 198022430),
            ('chr4', 191154276),
            ('chr5', 180915260),
            ('chr6', 171115067),
            ('chr7', 159138663),
            ('chr8', 146364022),
            ('chr9', 141213431),
            ('chr10', 135534747),
            ('chr11', 135006516),
            ('chr12', 133851895),
            ('chr13', 115169878),
            ('chr14', 107349540),
            ('chr15', 102531392),
            ('chr16', 90354753),
            ('chr17', 81195210),
            ('chr18', 78077248),
            ('chr19', 59128983),
            ('chr20', 63025520),
            ('chr21', 48129895),
            ('chr22', 51304566),
            ('chrX', 155270560)])

    @classmethod
    def is_named(cls, name):
        return name == "hg19"

def build(name):
    for cls in Assembly.__subclasses__():
        if cls.is_named(name):
            return cls()
    raise ValueError
