from anadama_workflows import pipelines

def eq(a, b):
    assert a == b

def deepeq(answers, result):
    assert (a==b for a, b in zip(answers, result))

def test__regex_filter():
    a = ['1.fastq', '2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0]], [a[1]], []), result

    a = ['r1.fastq', 'R2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0]], [a[1]], []), result

    a = ['sample_a.r1.fastq', 'sample_b.r2.fastq', 'sample_c.r1.fastq',
         'sample_b.r1.fastq', 'sample_c.r2.fastq', 'sample_a.r2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0], a[3], a[2]], [a[1], a[5], a[4]], []), result

    a = ['sample.a.r1.fastq', 'sample_b.R2.fastq', 'sample.c.r1.fastq',
         'sample-b.r1.fastq', 'sample_c.r2.fastq', 'sample-a.R2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0], a[3], a[2]], [a[1], a[5], a[4]], []), result


def test__to_merged():
    yield eq, pipelines._to_merged("koji_R1.fastq", "foo"), "koji_foo.fastq"
    yield eq, pipelines._to_merged("koji_R1.fastq"), "koji_merged.fastq"
    yield eq, pipelines._to_merged("kojiR1.fastq", "foo"), "kojifoo.fastq"
    yield eq, pipelines._to_merged("koji1.fastq"), "kojimerged.fastq"
    yield eq, pipelines._to_merged("koji1.fastq"), "kojimerged.fastq"
    yield eq, pipelines._to_merged("2.fastq", "stitched"), "stitched.fastq"
    
    
