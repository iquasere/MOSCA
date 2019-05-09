# -*- coding: utf-8 -*-
"""
MOSCA's MetaProteomics class for performing MetaProteomics Analysis

By João Sequeira

May 2019
"""

class UniprotMapping:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
        uniprot_columns = {'Entry':'id',                                        # From https://www.uniprot.org/help/uniprotkb_column_names
                    'Entry name':'entry name',
                    'Gene names':'genes',
                    'Gene names (primary)':'genes(PREFERRED)',
                    'Gene names (synonym)':'genes(ALTERNATIVE)',
                    'Gene names (ordered locus)':'genes(OLN)',
                    'Gene names (ORF)':'genes(ORF)',
                    'Organism':'organism',
                    'Organism ID':'organism-id',
                    'Protein names':'protein names',
                    'Proteomes':'proteome',
                    'Taxonomic lineage':'lineage(ALL)',
                    'Virus hosts':'virus hosts',
                    'Fragment':'fragment',
                    'Gene encoded by':'encodedon',
                    'Alternative products':'comment(ALTERNATIVE PRODUCTS)',
                    'Erroneous gene model prediction':'comment(ERRONEOUS GENE MODEL PREDICTION)',
                    'Erroneous initiation':'comment(ERRONEOUS INITIATION)',
                    'Erroneous termination':'comment(ERRONEOUS TERMINATION)',
                    'Erroneous translation':'comment(ERRONEOUS TRANSLATION)',
                    'Frameshift':'comment(FRAMESHIFT)',
                    'Mass spectrometry':'comment(MASS SPECTROMETRY)',
                    'Polymorphism':'comment(POLYMORPHISM)',
                    'RNA editing':'comment(RNA EDITING)',
                    'Sequence caution':'comment(SEQUENCE CAUTION)',
                    'Length':'length',
                    'Mass':'mass',
                    'Sequence':'sequence',
                    'Alternative sequence':'feature(ALTERNATIVE SEQUENCE)',
                    'Natural variant':'feature(NATURAL VARIANT)',
                    'Non-adjacent residues':'feature(NON ADJACENT RESIDUES)',
                    'Non-standard residue':'feature(NON STANDARD RESIDUE)',
                    'Non-terminal residue':'feature(NON TERMINAL RESIDUE)',
                    'Sequence conflict':'feature(SEQUENCE CONFLICT)',
                    'Sequence uncertainty':'feature(SEQUENCE UNCERTAINTY)',
                    'Sequence version':'version(sequence)',
                    'EC number':'ec',
                    'Absorption':'comment(ABSORPTION)',
                    'Catalytic activity':'comment(CATALYTIC ACTIVITY)',
                    'ChEBI':'chebi',
                    'ChEBI (Catalytic activity)':'chebi(Catalytic activity)',
                    'ChEBI (Cofactor)':'chebi(Cofactor)',
                    'ChEBI IDs':'chebi-id',
                    'Cofactor':'comment(COFACTOR)',
                    'Enzyme regulation':'comment(ENZYME REGULATION)',
                    'Function[CC]':'comment(FUNCTION)',
                    'Kinetics':'comment(KINETICS)',
                    'Pathway':'comment(PATHWAY)',
                    'Redox potential':'comment(REDOX POTENTIAL)',
                    'Temperature dependence':'comment(TEMPERATURE DEPENDENCE)',
                    'pH dependence':'comment(PH DEPENDENCE)',
                    'Active site':'feature(ACTIVE SITE)',
                    'Binding site':'feature(BINDING SITE)',
                    'DNA binding':'feature(DNA BINDING)',
                    'Metal binding':'feature(METAL BINDING)',
                    'Nucleotide binding':'feature(NP BIND)',
                    'Site':'feature(SITE)',
                    'Annotation score':'annotation score',
                    'Features':'features',
                    'Caution':'comment(CAUTION)',
                    'Miscellaneous[CC]':'comment(MISCELLANEOUS)',
                    'Keywords':'keywords',
                    'Matched text':'context',
                    'Protein existence':'existence',
                    'Tools':'tools',
                    'Reviewed':'reviewed',
                    'Subunit structure[CC]':'comment(SUBUNIT)',
                    'Interacts with':'interactor',
                    'Developmental stage':'comment(DEVELOPMENTAL STAGE)',
                    'Induction':'comment(INDUCTION)',
                    'Tissue specificity':'comment(TISSUE SPECIFICITY)',
                    'Gene ontology (GO)':'go',
                    'Gene ontology (biological process)':'go(biological process)',
                    'Gene ontology (molecular function)':'go(molecular function)',
                    'Gene ontology (cellular component)':'go(cellular component)',
                    'Gene ontology IDs':'go-id',
                    'Allergenic properties':'comment(ALLERGEN)',
                    'Biotechnological use':'comment(BIOTECHNOLOGY)',
                    'Disruption phenotype':'comment(DISRUPTION PHENOTYPE)',
                    'Involvement in disease':'comment(DISEASE)',
                    'Pharmaceutical use':'comment(PHARMACEUTICAL)',
                    'Toxic dose':'comment(TOXIC DOSE)',
                    'Subcellular location[CC]':'comment(SUBCELLULAR LOCATION)',
                    'Intramembrane':'feature(INTRAMEMBRANE)',
                    'Topological domain':'feature(TOPOLOGICAL DOMAIN)',
                    'Transmembrane':'feature(TRANSMEMBRANE)',
                    'Post-translational modification':'comment(PTM)',
                    'Chain':'feature(CHAIN)',
                    'Cross-link':'feature(CROSS LINK)',
                    'Disulfide bond':'feature(DISULFIDE BOND)',
                    'Glycosylation':'feature(GLYCOSYLATION)',
                    'Initiator methionine':'feature(INITIATOR METHIONINE)',
                    'Lipidation':'feature(LIPIDATION)',
                    'Modified residue':'feature(MODIFIED RESIDUE)',
                    'Peptide':'feature(PEPTIDE)',
                    'Propeptide':'feature(PROPEPTIDE)',
                    'Signal peptide':'feature(SIGNAL)',
                    'Transit peptide':'feature(TRANSIT)',
                    '3D':'3d',
                    'Beta strand':'feature(BETA STRAND)',
                    'Helix':'feature(HELIX)',
                    'Turn':'feature(TURN)',
                    'Mapped PubMed ID':'citationmapping',
                    'PubMed ID':'citation',
                    'Date of creation':'created',
                    'Date of last modification':'last-modified',
                    'Date of last sequence modification':'sequence-modified',
                    'Entry version':'version(entry)',
                    'Domain[CC]':'comment(DOMAIN)',
                    'Sequence similarities':'comment(SIMILARITY)',
                    'Protein families':'families',
                    'Coiled coil':'feature(COILED COIL)',
                    'Compositional bias':'feature(COMPOSITIONAL BIAS)',
                    'Domain[FT]':'feature(DOMAIN EXTENT)',
                    'Motif':'feature(MOTIF)',
                    'Region':'feature(REGION)',
                    'Repeat':'feature(REPEAT)',
                    'Zinc finger':'feature(ZINC FINGER)',
                    'Taxonomic lineage (all)':'lineage(all)',
                    'Taxonomic lineage (SUPERKINGDOM)':'lineage(SUPERKINGDOM)',
                    'Taxonomic lineage (KINGDOM)':'lineage(KINGDOM)',
                    'Taxonomic lineage (SUBKINGDOM)':'lineage(SUBKINGDOM)',
                    'Taxonomic lineage (SUPERPHYLUM)':'lineage(SUPERPHYLUM)',
                    'Taxonomic lineage (PHYLUM)':'lineage(PHYLUM)',
                    'Taxonomic lineage (SUBPHYLUM)':'lineage(SUBPHYLUM)',
                    'Taxonomic lineage (SUPERCLASS)':'lineage(SUPERCLASS)',
                    'Taxonomic lineage (CLASS)':'lineage(CLASS)',
                    'Taxonomic lineage (SUBCLASS)':'lineage(SUBCLASS)',
                    'Taxonomic lineage (INFRACLASS)':'lineage(INFRACLASS)',
                    'Taxonomic lineage (SUPERORDER)':'lineage(SUPERORDER)',
                    'Taxonomic lineage (ORDER)':'lineage(ORDER)',
                    'Taxonomic lineage (SUBORDER)':'lineage(SUBORDER)',
                    'Taxonomic lineage (INFRAORDER)':'lineage(INFRAORDER)',
                    'Taxonomic lineage (PARVORDER)':'lineage(PARVORDER)',
                    'Taxonomic lineage (SUPERFAMILY)':'lineage(SUPERFAMILY)',
                    'Taxonomic lineage (FAMILY)':'lineage(FAMILY)',
                    'Taxonomic lineage (SUBFAMILY)':'lineage(SUBFAMILY)',
                    'Taxonomic lineage (TRIBE)':'lineage(TRIBE)',
                    'Taxonomic lineage (SUBTRIBE)':'lineage(SUBTRIBE)',
                    'Taxonomic lineage (GENUS)':'lineage(GENUS)',
                    'Taxonomic lineage (SUBGENUS)':'lineage(SUBGENUS)',
                    'Taxonomic lineage (SPECIES GROUP)':'lineage(SPECIES GROUP)',
                    'Taxonomic lineage (SPECIES SUBGROUP)':'lineage(SPECIES SUBGROUP)',
                    'Taxonomic lineage (SPECIES)':'lineage(SPECIES)',
                    'Taxonomic lineage (SUBSPECIES)':'lineage(SUBSPECIES)',
                    'Taxonomic lineage (VARIETAS)':'lineage(VARIETAS)',
                    'Taxonomic lineage (FORMA)':'lineage(FORMA)',
                    'Taxonomic identifier (all)':'lineage-id(all)',
                    'Taxonomic identifier (SUPERKINGDOM)':'lineage-id(SUPERKINGDOM)',
                    'Taxonomic identifier (KINGDOM)':'lineage-id(KINGDOM)',
                    'Taxonomic identifier (SUBKINGDOM)':'lineage-id(SUBKINGDOM)',
                    'Taxonomic identifier (SUPERPHYLUM)':'lineage-id(SUPERPHYLUM)',
                    'Taxonomic identifier (PHYLUM)':'lineage-id(PHYLUM)',
                    'Taxonomic identifier (SUBPHYLUM)':'lineage-id(SUBPHYLUM)',
                    'Taxonomic identifier (SUPERCLASS)':'lineage-id(SUPERCLASS)',
                    'Taxonomic identifier (CLASS)':'lineage-id(CLASS)',
                    'Taxonomic identifier (SUBCLASS)':'lineage-id(SUBCLASS)',
                    'Taxonomic identifier (INFRACLASS)':'lineage-id(INFRACLASS)',
                    'Taxonomic identifier (SUPERORDER)':'lineage-id(SUPERORDER)',
                    'Taxonomic identifier (ORDER)':'lineage-id(ORDER)',
                    'Taxonomic identifier (SUBORDER)':'lineage-id(SUBORDER)',
                    'Taxonomic identifier (INFRAORDER)':'lineage-id(INFRAORDER)',
                    'Taxonomic identifier (PARVORDER)':'lineage-id(PARVORDER)',
                    'Taxonomic identifier (SUPERFAMILY)':'lineage-id(SUPERFAMILY)',
                    'Taxonomic identifier (FAMILY)':'lineage-id(FAMILY)',
                    'Taxonomic identifier (SUBFAMILY)':'lineage-id(SUBFAMILY)',
                    'Taxonomic identifier (TRIBE)':'lineage-id(TRIBE)',
                    'Taxonomic identifier (SUBTRIBE)':'lineage-id(SUBTRIBE)',
                    'Taxonomic identifier (GENUS)':'lineage-id(GENUS)',
                    'Taxonomic identifier (SUBGENUS)':'lineage-id(SUBGENUS)',
                    'Taxonomic identifier (SPECIES GROUP)':'lineage-id(SPECIES GROUP)',
                    'Taxonomic identifier (SPECIES SUBGROUP)':'lineage-id(SPECIES SUBGROUP)',
                    'Taxonomic identifier (SPECIES)':'lineage-id(SPECIES)',
                    'Taxonomic identifier (SUBSPECIES)':'lineage-id(SUBSPECIES)',
                    'Taxonomic identifier (VARIETAS)':'lineage-id(VARIETAS)',
                    'Taxonomic identifier (FORMA)':'lineage-id(FORMA)',
                    'db_abbrev':'database(db_abbrev)'}
        
        uniprot_databases = {'Allergome; a platform for allergen knowledge':'Allergome',
                    'ArachnoServer: Spider toxin database':'ArachnoServer',
                    'Arabidopsis Information Portal':'Araport',
                    'Bgee dataBase for Gene Expression Evolution':'Bgee',
                    'BindingDB database of measured binding affinities':'BindingDB',
                    'BioCyc Collection of Pathway/Genome Databases':'BioCyc',
                    'The Biological General Repository for Interaction Datasets (BioGrid)':'BioGrid',
                    'BioMuta curated single-nucleotide variation and disease association database':'BioMuta',
                    'BRENDA Comprehensive Enzyme Information System':'BRENDA',
                    'CarbonylDB database of protein carbonylation sites':'CarbonylDB',
                    'Carbohydrate-Active enZymes':'CAZy',
                    'The Consensus CDS (CCDS) project':'CCDS',
                    'Conserved Domains Database':'CDD',
                    'Candida Genome Database':'CGD',
                    'ChEMBL database of bioactive drug-like small molecules':'ChEMBL',
                    'ChiTaRS: a database of human, mouse and fruit fly chimeric transcripts and RNA-sequencing data':'ChiTaRS',
                    'CollecTF database of bacterial transcription factor binding sites':'CollecTF',
                    'ComplexPortal: manually curated resource of macromolecular complexes':'ComplexPortal',
                    '2-DE database at Universidad Complutense de Madrid':'COMPLUYEAST-2DPAGE',
                    'ConoServer: Cone snail toxin database':'ConoServer',
                    'CORUM comprehensive resource of mammalian protein complexes':'CORUM',
                    'Comparative Toxicogenomics Database':'CTD',
                    'Database of single nucleotide polymorphism':'dbSNP',
                    'DNA Data Bank of Japan; a nucleotide sequence database':'DDBJ',
                    'DEPOD human dephosphorylation database':'DEPOD',
                    'Dictyostelium discoideum online informatics resource':'dictyBase',
                    'Database of interacting proteins':'DIP',
                    'DisGeNET':'DisGeNET',
                    'Database of protein disorder':'DisProt',
                    'Domain mapping of disease mutations (DMDM)':'DMDM',
                    'The DNASU plasmid repository':'DNASU',
                    'DOSAC-COBS 2D-PAGE database':'DOSAC-COBS-2DPAGE',
                    'Drug and drug target database':'DrugBank',
                    'EchoBASE - an integrated post-genomic database for E. coli':'EchoBASE',
                    'Escherichia coli strain K12 genome database':'EcoGene',
                    'evolutionary genealogy of genes: Non-supervised Orthologous Groups':'eggNOG',
                    'The Eukaryotic Linear Motif resource for Functional Sites in Proteins':'ELM',
                    'EMBL nucleotide sequence database':'EMBL',
                    'Ensembl eukaryotic genome annotation project':'Ensembl',
                    'Ensembl bacterial and archaeal genome annotation project':'EnsemblBacteria',
                    'Ensembl fungal genome annotation project':'EnsemblFungi',
                    'Ensembl metazoan genome annotation project':'EnsemblMetazoa',
                    'Ensembl plant genome annotation project':'EnsemblPlants',
                    'Ensembl protists genome annotation project':'EnsemblProtists',
                    'Enzyme nomenclature database':'ENZYME',
                    'Encyclopedia of Proteome Dynamics':'EPD',
                    'ESTHER database of the Alpha/Beta-hydrolase fold superfamily of proteins':'ESTHER',
                    'European Hepatitis C Virus Database':'euHCVdb',
                    'Eukaryotic Pathogen Database Resources':'EuPathDB',
                    'Relative evolutionary importance of amino acids within a protein sequence':'EvolutionaryTrace',
                    'ExpressionAtlas, Differential and Baseline Expression':'ExpressionAtlas',
                    'Drosophila genome database':'FlyBase',
                    'GenAtlas: human gene database':'GenAtlas',
                    'GenBank nucleotide sequence database':'GenBank',
                    'Gene3D Structural and Functional Annotation of Protein Families':'Gene3D',
                    'GeneCards: human genes, protein and diseases':'GeneCards',
                    'GeneDB pathogen genome database from Sanger Institute':'GeneDB',
                    'Database of genes from NCBI RefSeq genomes':'GeneID',
                    'GeneReviews a resource of expert-authored, peer-reviewed disease descriptions.':'GeneReviews',
                    'Ensembl GeneTree':'GeneTree',
                    'Genevisible search portal to normalized and curated expression data from Genevestigator':'Genevisible',
                    'The Gene Wiki collection of pages on human genes and proteins':'GeneWiki',
                    'Database of phenotypes from RNA interference screens in Drosophila and Homo sapiens':'GenomeRNAi',
                    'GlyConnect protein glycosylation platform':'GlyConnect',
                    'Gene Ontology':'GO',
                    'Information system for G protein-coupled receptors (GPCRs)':'GPCRDB',
                    'Gramene; a comparative resource for plants':'Gramene',
                    'IUPHAR/BPS Guide to PHARMACOLOGY':'GuidetoPHARMACOLOGY',
                    'H-Invitational Database, human transcriptome db':'H-InvDB',
                    'HAMAP database of protein families':'HAMAP',
                    'Human Gene Nomenclature Database':'HGNC',
                    'The HOGENOM Database of Homologous Genes from Fully Sequenced Organisms':'HOGENOM',
                    'Human Protein Atlas':'HPA',
                    'Human Unidentified Gene-Encoded large proteins database':'HUGE',
                    'The international ImMunoGeneTics information system':'IMGT_GENE-DB',
                    'InParanoid: Eukaryotic Ortholog Groups':'InParanoid',
                    'Protein interaction database and analysis system':'IntAct',
                    'Integrated resource of protein families, domains and functional sites':'InterPro',
                    'iPTMnet integrated resource for PTMs in systems biology context':'iPTMnet',
                    'jPOST - Japan Proteome Standard Repository/Database':'jPOST',
                    'KEGG: Kyoto Encyclopedia of Genes and Genomes':'KEGG',
                    'KEGG Orthology (KO)':'KO',
                    'Legionella pneumophila genome database':'LegioList',
                    'Mycobacterium leprae genome database':'Leproma',
                    'Maize Genetics and Genomics Database':'MaizeGDB',
                    'MalaCards human disease database':'MalaCards',
                    'MaxQB - The MaxQuant DataBase':'MaxQB',
                    'MEROPS protease database':'MEROPS',
                    'Mouse genome database (MGD) from Mouse Genome Informatics (MGI)':'MGI',
                    'Microbial advanced database':'Micado',
                    'Online Mendelian Inheritance in Man (OMIM)':'MIM',
                    'Molecular INTeraction database':'MINT',
                    'MobiDB: a database of protein disorder and mobility annotations':'MobiDB',
                    'Database of comparative protein structure models':'ModBase',
                    'MoonDB Database of extreme multifunctional and moonlighting proteins':'MoonDB',
                    'MoonProt database of moonlighting proteins':'MoonProt',
                    'mycoCLAP, a database of fungal genes encoding lignocellulose-active proteins':'mycoCLAP',
                    'neXtProt; the human protein knowledge platform':'neXtProt',
                    'USC-OGP 2-DE database':'OGP',
                    'Identification of Orthologs from Complete Genome Data':'OMA',
                    'Open Targets':'OpenTargets',
                    'Orphanet; a database dedicated to information on rare diseases and orphan drugs':'Orphanet',
                    'Database of Orthologous Groups':'OrthoDB',
                    'The PANTHER Classification System':'PANTHER',
                    'Pathosystems Resource Integration Center (PATRIC)':'PATRIC',
                    'PaxDb, a database of protein abundance averages across all three domains of life':'PaxDb',
                    'Protein Data Bank Europe':'PDB',
                    'Protein Data Bank Japan':'PDBj',
                    'PDBsum; at-a-glance overview of macromolecular structures':'PDBsum',
                    'PeptideAtlas':'PeptideAtlas',
                    'PeroxiBase, a peroxidase database':'PeroxiBase',
                    'Pfam protein domain database':'Pfam',
                    'The Pharmacogenetics and Pharmacogenomics Knowledge Base':'PharmGKB',
                    'Comprehensive resource for the study of protein post-translational modifications (PTMs) in human, mouse and rat.':'PhosphoSitePlus',
                    'Database for complete collections of gene phylogenies':'PhylomeDB',
                    'Protein sequence database of the Protein Information Resource':'PIR',
                    'PIRSF; a whole-protein classification database':'PIRSF',
                    'CutDB - Proteolytic event database':'PMAP-CutDB',
                    'Schizosaccharomyces pombe database':'PomBase',
                    'PRoteomics IDEntifications database':'PRIDE',
                    'Protein Motif fingerprint database; a protein domain database':'PRINTS',
                    'Protein Ontology':'PRO',
                    'ProDom; a protein domain database':'ProDom',
                    'Protein Mass spectra EXtraction':'ProMEX',
                    'PROSITE; a protein domain and family database':'PROSITE',
                    'Proteomes':'Proteomes',
                    'ProteomicsDB human proteome resource':'ProteomicsDB',
                    'ProtoNet; Automatic hierarchical classification of proteins':'ProtoNet',
                    'Pseudomonas genome database':'PseudoCAP',
                    'Protein Data Bank RCSB':'RCSB-PDB',
                    'Reactome - a knowledgebase of biological pathways and processes':'Reactome',
                    'Restriction enzymes and methylases database':'REBASE',
                    'NCBI Reference Sequences':'RefSeq',
                    'REPRODUCTION-2DPAGE':'REPRODUCTION-2DPAGE',
                    'Rat genome database':'RGD',
                    'Rodent Unidentified Gene-Encoded large proteins database':'Rouge',
                    'SABIO-RK: Biochemical Reaction Kinetics Database':'SABIO-RK',
                    'The Structural Biology Knowledgebase':'SBKB',
                    'Structure-Function Linkage Database':'SFLD',
                    'Saccharomyces Genome Database':'SGD',
                    'SignaLink: a signaling pathway resource with multi-layered regulatory networks':'SignaLink',
                    'SIGNOR Signaling Network Open Resource':'SIGNOR',
                    'Simple Modular Architecture Research Tool; a protein domain database':'SMART',
                    'SWISS-MODEL Repository - a database of annotated 3D protein structure models':'SMR',
                    'The Stanford Online Universal Resource for Clones and ESTs':'SOURCE',
                    'STRING: functional protein association networks':'STRING',
                    'Superfamily database of structural and functional annotation':'SUPFAM',
                    'Two-dimensional polyacrylamide gel electrophoresis database from the Geneva University Hospital':'SWISS-2DPAGE',
                    'SWISS-MODEL Interactive Workspace':'SWISS-MODEL-Workspace',
                    'SwissLipids knowledge resource for lipid biology':'SwissLipids',
                    'SwissPalm database of S-palmitoylation events':'SwissPalm',
                    'The Arabidopsis Information Resource':'TAIR',
                    'Transport Classification Database':'TCDB',
                    'TIGRFAMs; a protein family database':'TIGRFAMs',
                    'Consortium for Top Down Proteomics':'TopDownProteomics',
                    'TreeFam database of animal gene trees':'TreeFam',
                    'Mycobacterium tuberculosis strain H37Rv genome database':'TubercuList',
                    'University College Dublin 2-DE Proteome Database':'UCD-2DPAGE',
                    'UCSC genome browser':'UCSC',
                    'UniCarbKB; an annotated and curated database of glycan structures':'UniCarbKB',
                    'UniLectin database of carbohydrate-binding proteins':'UniLectin',
                    'UniPathway: a resource for the exploration and annotation of metabolic pathways':'UniPathway',
                    'Bioinformatics Resource for Invertebrate Vectors of Human Pathogens':'VectorBase',
                    'Vertebrate Gene Nomenclature Database':'VGNC',
                    'WormBase ParaSite':'WBParaSite',
                    'The World-2DPAGE database':'World-2DPAGE',
                    'WormBase':'WormBase',
                    'Xenopus laevis and tropicalis biology and genomics resource':'Xenbase',
                    'Zebrafish Information Network genome database':'ZFIN'}
        
        def string4mapping(self, columns = None, databases = None):
            if columns is None and databases is None:                           # Sets to defaults, it's dirty but best way I found
                columns = ['Entry', 'Gene names', 'Protein names',
                'EC number', 'Function[CC]', 'Pathway', 'Keywords', 'Protein existence', 
                'Gene ontology (GO)', 'Protein families', 'Taxonomic lineage (SUPERKINGDOM)', 
                'Taxonomic lineage (PHYLUM)', 'Taxonomic lineage (CLASS)', 
                'Taxonomic lineage (ORDER)', 'Taxonomic lineage (FAMILY)', 
                'Taxonomic lineage (GENUS)', 'Taxonomic lineage (SPECIES)'],
                databases = ['BioCyc Collection of Pathway/Genome Databases',
                 'BRENDA Comprehensive Enzyme Information System',
                 'Conserved Domains Database', 
                 'evolutionary genealogy of genes: Non-supervised Orthologous Groups',
                 'Ensembl eukaryotic genome annotation project',
                 'Integrated resource of protein families, domains and functional sites',
                 'KEGG: Kyoto Encyclopedia of Genes and Genomes',
                 'KEGG Orthology (KO)', 'Pfam protein domain database',
                 'Reactome - a knowledgebase of biological pathways and processes',
                 'NCBI Reference Sequences',
                 'UniPathway: a resource for the exploration and annotation of metabolic pathways']
            result = ','.join(columns)
            if len(databases) > 0:
                result += 'database(' + '),database('.join([uniprot_databases(db) for db in databases]) + ')'
            return result