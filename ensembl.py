import wrenlab

class ensembl_loci(object):
    def __init__(self,data=None):
        self.data = data
        probes = "GPL13534-11288.txt"
        probes = pd.read_csv(probes,skiprows=37,sep="\t",index_col=0)
        probes = probes[["Name","Chromosome_36","RANGE_START","RANGE_END"]]
        probes.columns = ["ProbeID","chrom","start","end"]
        probes["strand"] = ""
        probes = to_bedtool(probes)
        self.probes = probes
        loci = wrenlab.ensembl.loci(9606,version=75)
        #save df version to change later if you need it 
        self.loci = loci 
        self.genes = to_bedtool(loci.reset_index())
        #df = gene.ensembl(9606)
        df = "ensembl_9606.csv"
        df = pd.read_csv(df,sep=",",header=0,index_col=0)
        self.ensembl_entrez = dict(zip(df["Ensembl Gene ID"], df["GeneID"]))
        
    #window relative to the start (l) to end of gene (r)
    def get_gene_window(self,l,r):
        intersection = self.genes.window(self.probes, l=l, r=r, sw=True)
        locus_gene = {obj.fields[8]:self.ensembl_entrez[obj.fields[3]] for obj in intersection if obj.fields[3] in self.ensembl_entrez }
        return(locus_gene)
    
    def gene_window_collapse(self,l,r,method="mean"):
        locus_gene = self.get_gene_window(l=l,r=r)
        out = self.data.groupby(locus_gene).apply(lambda x : collapse_genes(x,method=method))
        return(out)
    
    def get_tss_window(self,l,r):
        genes = self.loci
        genes["end"] = genes["start"]
        genes = to_bedtool(genes.reset_index())
        intersection = genes.window(self.probes, l=l, r=r, sw=True)
        locus_gene = {obj.fields[8]:self.ensembl_entrez[obj.fields[3]] for obj in intersection if obj.fields[3] in self.ensembl_entrez }
        #locus_gene = self.get_gene_window(l,r,genes)
        return(locus_gene)
    
    def tss_window_collapse(self,l,r,method="mean"):
        locus_gene = self.get_tss_window(l,r)
        out = self.data.groupby(locus_gene).apply(lambda x : collapse_genes(x,method=method))
        return(out)

    def tss20002000(self,method="mean"):
        out = self.tss_window_collapse(l=2000,r=2000,method=method)
        return(out)
    
    def tss1500(self,method="mean"):
        out = self.tss_window_collapse(l=1500,r=0,method=method)
        return(out)
    
    def tss200(self,method="mean"):
        out = self.tss_window_collapse(l=200,r=0,method=method)
        return(out)
    
    def tss5050(self,method="mean"):
        out = self.tss_window_collapse(l=50,r=50,method=method)
        return(out)
    
    def gb(self,method="mean"):
        out = self.gene_window_collapse(l=0,r=0,method=method)
        return(out)
    
    def extended_gb(self,method="mean"):
        out = self.gene_window_collapse(l=1500,r=500,method=method)
        return(out)
    
    def exon():
        pass
    
    def intron():
        pass    
