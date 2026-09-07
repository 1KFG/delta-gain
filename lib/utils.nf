// Shared Groovy helpers for DeltaGain processes.
//
// Adapted from (not copied verbatim from) Fungi_BFD/nextflow/modules/common/utils.nf:
// that file's skani*For()/aniTimeFor() helpers scale for per-species/genus groups
// capped around 500 genomes, and aniTimeFor() there also carries a k8s-pilot-only
// 6h hard cap that would be wrong here. DeltaGain's genome-classification stage runs
// ONE global group (~22k BFD genomes + reference genomes, tens of thousands), so the
// tiers below are sized for that instead. No k8s path exists for this project yet,
// so the capCpus/capMemGB k8s-clamp indirection is dropped too.

// mash sketch/dist scale roughly linearly in genome count; give it real headroom
// at the tens-of-thousands scale this pipeline actually runs at.
def mashCpusFor(int n) {
    n > 15000 ? 64 : n > 5000 ? 32 : n > 500 ? 16 : 8
}

def mashMemoryFor(int n, int attempt = 1) {
    def baseGB  = n > 15000 ? 256 : n > 5000 ? 128 : n > 500 ? 48 : 16
    def floorGB = attempt >= 3 ? 384 : attempt == 2 ? 192 : baseGB
    Math.max(baseGB, floorGB).toString() + ' GB'
}

def mashTimeFor(int n, int attempt = 1) {
    def baseH  = n > 15000 ? 12 : n > 5000 ? 6 : n > 500 ? 2 : 1
    def floorH = attempt >= 2 ? 24 : baseH
    Math.max(baseH, floorH).toString() + ' h'
}

// Per-component skani triangle refinement: components are expected to be much
// smaller than the global mash group (a component is one candidate species/strain
// cluster), but a few large well-sampled species (e.g. Fusarium, Aspergillus
// species complexes) can still run into the thousands — same shape as the
// reference project's skaniMemoryFor/aniTimeFor, just inlined here so this
// project doesn't depend on a sibling project's file layout.
def skaniCpusFor(int n) {
    n > 1000 ? 32 : n > 200 ? 16 : 8
}

def skaniMemoryFor(int n, int attempt = 1) {
    def baseGB  = n > 1000 ? 96 : n > 200 ? 32 : 16
    def floorGB = attempt >= 3 ? 192 : attempt == 2 ? 48 : baseGB
    Math.max(baseGB, floorGB).toString() + ' GB'
}

def skaniTimeFor(int n, int attempt = 1) {
    def baseH  = n > 1000 ? 8 : n > 200 ? 3 : 1
    def floorH = attempt >= 2 ? 24 : baseH
    Math.max(baseH, floorH).toString() + ' h'
}
