#!/usr/bin/env python

from Bio.KEGG.REST import kegg_get, kegg_link, kegg_list
from Bio.KEGG.KGML import KGML_parser, KGML_pathway
from Bio.Graphics.KGML_vis import KGMLCanvas
from matplotlib import cm
from matplotlib.colors import to_hex
from io import StringIO
import matplotlib.pyplot as plt
import re, pandas as pd, numpy as np, os
from progressbar import ProgressBar

__author__ = "Tiago Oliveira"
__credits__ = ["Tiago Oliveira", "Joao Sequeira"]
__version__ = "1.0"
__maintainer__ = "Joao Sequeira"
__email__ = "maildosequeira@gmail.com"
__status__ = "Production"
################################################################################
class KeggPathway:
    '''
    This class concerns with the storage of information manually retrieved from 
    KEGG documentation for reference, and the conversion of certain IDs to 
    different IDs from the KEGG DB
    '''
    
    def __init__(self, **kwargs):
        '''
        Initialize object
        :param **kwargs - none is necessary for this class
        '''
        
        self.__dict__ = kwargs
        
        self.maps = {
            # 1. Metabolism
            # 1.0 Global and overview maps
            '01100':'Metabolic pathways',
            '01110':'Biosynthesis of secondary metabolites',
            '01120':'Microbial metabolism in diverse environments',
            '01130':'Biosynthesis of antibiotics',
            '01200':'Carbon metabolism',
            '01210':'2-Oxocarboxylic acid metabolism',
            '01212':'Fatty acid metabolism',
            '01230':'Biosynthesis of amino acids',
            '01220':'Degradation of aromatic compounds',
            # 1.1 Carbohydrate metabolism
            '00010':'Glycolysis / Gluconeogenesis',
            '00020':'Citrate cycle (TCA cycle)',
            '00030':'Pentose phosphate pathway',
            '00040':'Pentose and glucuronate interconversions',
            '00051':'Fructose and mannose metabolism',
            '00052':'Galactose metabolism',
            '00053':'Ascorbate and aldarate metabolism',
            '00500':'Starch and sucrose metabolism',
            '00520':'Amino sugar and nucleotide sugar metabolism',
            '00620':'Pyruvate metabolism',
            '00630':'Glyoxylate and dicarboxylate metabolism',
            '00640':'Propanoate metabolism',
            '00650':'Butanoate metabolism',
            '00660':'C5-Branched dibasic acid metabolism',
            '00562':'Inositol phosphate metabolism',
            # 1.2 Energy metabolism
            '00190':'Oxidative phosphorylation',
            '00195':'Photosynthesis',
            '00196':'Photosynthesis - antenna proteins',
            '00710':'Carbon fixation in photosynthetic organisms',
            '00720':'Carbon fixation pathways in prokaryotes',
            '00680':'Methane metabolism',
            '00910':'Nitrogen metabolism',
            '00920':'Sulfur metabolism',
            # 1.3 Lipid metabolism
            '00061':'Fatty acid biosynthesis',
            '00062':'Fatty acid elongation',
            '00071':'Fatty acid degradation',
            '00072':'Synthesis and degradation of ketone bodies',
            '00073':'Cutin, suberine and wax biosynthesis',
            '00100':'Steroid biosynthesis',
            '00120':'Primary bile acid biosynthesis',
            '00121':'Secondary bile acid biosynthesis',
            '00140':'Steroid hormone biosynthesis',
            '00561':'Glycerolipid metabolism',
            '00564':'Glycerophospholipid metabolism',
            '00565':'Ether lipid metabolism',
            '00600':'Sphingolipid metabolism',
            '00590':'Arachidonic acid metabolism',
            '00591':'Linoleic acid metabolism',
            '00592':'alpha-Linolenic acid metabolism',
            '01040':'Biosynthesis of unsaturated fatty acids',
            # 1.4 Nucleotide metabolism
            '00230':'Purine metabolism',
            '00240':'Pyrimidine metabolism',
            # 1.5 Amino acid metabolism
            '00250':'Alanine, aspartate and glutamate metabolism',
            '00260':'Glycine, serine and threonine metabolism',
            '00270':'Cysteine and methionine metabolism',
            '00280':'Valine, leucine and isoleucine degradation',
            '00290':'Valine, leucine and isoleucine biosynthesis',
            '00300':'Lysine biosynthesis','00310':'Lysine degradation',
            '00220':'Arginine biosynthesis',
            '00330':'Arginine and proline metabolism',
            '00340':'Histidine metabolism',
            '00350':'Tyrosine metabolism',
            '00360':'Phenylalanine metabolism',
            '00380':'Tryptophan metabolism',
            '00400':'Phenylalanine, tyrosine and tryptophan biosynthesis',
            # 1.6 Metabolism of other amino acids
            '00410':'beta-Alanine metabolism',
            '00430':'Taurine and hypotaurine metabolism',
            '00440':'Phosphonate and phosphinate metabolism',
            '00450':'Selenocompound metabolism',
            '00460':'Cyanoamino acid metabolism',
            '00471':'D-Glutamine and D-glutamate metabolism',
            '00472':'D-Arginine and D-ornithine metabolism',
            '00473':'D-Alanine metabolism',
            '00480':'Glutathione metabolism',
            # 1.7 Glycan biosynthesis and metabolism
            '00510':'N-Glycan biosynthesis',
            '00513':'Various types of N-glycan biosynthesis',
            '00512':'Mucin type O-glycan biosynthesis',
            '00515':'Mannose type O-glycan biosynthesis',
            '00514':'Other types of O-glycan biosynthesis',
            '00532':'Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate',
            '00534':'Glycosaminoglycan biosynthesis - heparan sulfate / heparin',
            '00533':'Glycosaminoglycan biosynthesis - keratan sulfate',
            '00531':'Glycosaminoglycan degradation',
            '00563':'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis',
            '00601':'Glycosphingolipid biosynthesis - lacto and neolacto series',
            '00603':'Glycosphingolipid biosynthesis - globo and isoglobo series',
            '00604':'Glycosphingolipid biosynthesis - ganglio series',
            '00540':'Lipopolysaccharide biosynthesis',
            '00550':'Peptidoglycan biosynthesis',
            '00511':'Other glycan degradation',
            '00571':'Lipoarabinomannan (LAM) biosynthesis',
            '00572':'Arabinogalactan biosynthesis - Mycobacterium',
            # 1.8 Metabolism of cofactors and vitamins
            '00730':'Thiamine metabolism',
            '00740':'Riboflavin metabolism',
            '00750':'Vitamin B6 metabolism',
            '00760':'Nicotinate and nicotinamide metabolism',
            '00770':'Pantothenate and CoA biosynthesis',
            '00780':'Biotin metabolism',
            '00785':'Lipoic acid metabolism',
            '00790':'Folate biosynthesis',
            '00670':'One carbon pool by folate',
            '00830':'Retinol metabolism',
            '00860':'Porphyrin and chlorophyll metabolism',
            '00130':'Ubiquinone and other terpenoid-quinone biosynthesis',
            # 1.9 Metabolism of terpenoids and polyketides
            '00900':'Terpenoid backbone biosynthesis',
            '00902':'Monoterpenoid biosynthesis',
            '00909':'Sesquiterpenoid and triterpenoid biosynthesis',
            '00904':'Diterpenoid biosynthesis',
            '00906':'Carotenoid biosynthesis',
            '00905':'Brassinosteroid biosynthesis',
            '00981':'Insect hormone biosynthesis',
            '00908':'Zeatin biosynthesis',
            '00903':'Limonene and pinene degradation',
            '00281':'Geraniol degradation',
            '01052':'Type I polyketide structures',
            '00522':'Biosynthesis of 12-, 14- and 16-membered macrolides',
            '01051':'Biosynthesis of ansamycins',
            '01059':'Biosynthesis of enediyne antibiotics',
            '01056':'Biosynthesis of type II polyketide backbone',
            '01057':'Biosynthesis of type II polyketide products',
            '00253':'Tetracycline biosynthesis',
            '00523':'Polyketide sugar unit biosynthesis',
            '01054':'Nonribosomal peptide structures',
            '01053':'Biosynthesis of siderophore group nonribosomal peptides',
            '01055':'Biosynthesis of vancomycin group antibiotics',
            # 1.10 Biosynthesis of other secondary metabolites
            '00940':'Phenylpropanoid biosynthesis',
            '00945':'Stilbenoid, diarylheptanoid and gingerol biosynthesis',
            '00941':'Flavonoid biosynthesis',
            '00944':'Flavone and flavonol biosynthesis',
            '00942':'Anthocyanin biosynthesis',
            '00943':'Isoflavonoid biosynthesis',
            '00901':'Indole alkaloid biosynthesis',
            '00403':'Indole diterpene alkaloid biosynthesis',
            '00950':'Isoquinoline alkaloid biosynthesis',
            '00960':'Tropane, piperidine and pyridine alkaloid biosynthesis',
            '01058':'Acridone alkaloid biosynthesis',
            '00232':'Caffeine metabolism',
            '00965':'Betalain biosynthesis',
            '00966':'Glucosinolate biosynthesis',
            '00402':'Benzoxazinoid biosynthesis',
            '00311':'Penicillin and cephalosporin biosynthesis',
            '00332':'Carbapenem biosynthesis',
            '00261':'Monobactam biosynthesis',
            '00331':'Clavulanic acid biosynthesis',
            '00521':'Streptomycin biosynthesis',
            '00524':'Neomycin, kanamycin and gentamicin biosynthesis',
            '00525':'Acarbose and validamycin biosynthesis',
            '00231':'Puromycin biosynthesis',
            '00401':'Novobiocin biosynthesis',
            '00404':'Staurosporine biosynthesis',
            '00405':'Phenazine biosynthesis',
            '00333':'Prodigiosin biosynthesis',
            '00254':'Aflatoxin biosynthesis',
            '00998':'Biosynthesis of secondary metabolites - other antibiotics New!',
            '00999':'Biosynthesis of secondary metabolites - unclassified',
            # 1.11 Xenobiotics biodegradation and metabolism
            '00362':'Benzoate degradation',
            '00627':'Aminobenzoate degradation',
            '00364':'Fluorobenzoate degradation',
            '00625':'Chloroalkane and chloroalkene degradation',
            '00361':'Chlorocyclohexane and chlorobenzene degradation',
            '00623':'Toluene degradation',
            '00622':'Xylene degradation',
            '00633':'Nitrotoluene degradation',
            '00642':'Ethylbenzene degradation',
            '00643':'Styrene degradation',
            '00791':'Atrazine degradation',
            '00930':'Caprolactam degradation',
            '00363':'Bisphenol degradation',
            '00621':'Dioxin degradation',
            '00626':'Naphthalene degradation',
            '00624':'Polycyclic aromatic hydrocarbon degradation',
            '00365':'Furfural degradation',
            '00984':'Steroid degradation',
            '00980':'Metabolism of xenobiotics by cytochrome P45',
            '00982':'2Drug metabolism - cytochrome P45',
            '00983':'3Drug metabolism - other enzymes',
            # 1.12 Chemical structure transformation maps
            '01010':'Overview of biosynthetic pathways',
            '01060':'Biosynthesis of plant secondary metabolites',
            '01061':'Biosynthesis of phenylpropanoids',
            '01062':'Biosynthesis of terpenoids and steroids',
            '01063':'Biosynthesis of alkaloids derived from shikimate pathway',
            '01064':'Biosynthesis of alkaloids derived from ornithine, lysine and nicotinic acid',
            '01065':'Biosynthesis of alkaloids derived from histidine and purine',
            '01066':'Biosynthesis of alkaloids derived from terpenoid and polyketide',
            '01070':'Biosynthesis of plant hormones',
            # 2. Genetic Information Processing
            # 2.1 Transcription
            '03020':'RNA polymerase',
            '03022':'Basal transcription factors',
            '03040':'Spliceosome',
            # 2.2 Translation
            '03010':'Ribosome',
            '00970':'Aminoacyl-tRNA biosynthesis',
            '03013':'RNA transport',
            '03015':'mRNA surveillance pathway',
            '03008':'Ribosome biogenesis in eukaryotes',
            # 2.3 Folding, sorting and degradation
            '03060':'Protein export',
            '04141':'Protein processing in endoplasmic reticulum',
            '04130':'SNARE interactions in vesicular transport',
            '04120':'Ubiquitin mediated proteolysis',
            '04122':'Sulfur relay system',
            '03050':'Proteasome',
            '03018':'RNA degradation',
            # 2.4 Replication and repair
            '03030':'DNA replication',
            '03410':'Base excision repair',
            '03420':'Nucleotide excision repair',
            '03430':'Mismatch repair',
            '03440':'Homologous recombination',
            '03450':'Non-homologous end-joining',
            '03460':'Fanconi anemia pathway',
            # 3. Environmental Information Processing
            # 3.1 Membrane transport
            '02010':'ABC transporters',
            '02060':'Phosphotransferase system (PTS)',
            '03070':'Bacterial secretion system',
            # 3.2 Signal transduction
            '02020':'Two-component system',
            '04014':'Ras signaling pathway',
            '04015':'Rap1 signaling pathway',
            '04010':'MAPK signaling pathway',
            '04013':'MAPK signaling pathway - fly',
            '04016':'MAPK signaling pathway - plant',
            '04011':'MAPK signaling pathway - yeast',
            '04012':'ErbB signaling pathway',
            '04310':'Wnt signaling pathway',
            '04330':'Notch signaling pathway',
            '04340':'Hedgehog signaling pathway',
            '04341':'Hedgehog signaling pathway - fly',
            '04350':'TGF-beta signaling pathway',
            '04390':'Hippo signaling pathway',
            '04391':'Hippo signaling pathway - fly',
            '04392':'Hippo signaling pathway - multiple species',
            '04370':'VEGF signaling pathway',
            '04371':'Apelin signaling pathway',
            '04630':'Jak-STAT signaling pathway',
            '04064':'NF-kappa B signaling pathway',
            '04668':'TNF signaling pathway',
            '04066':'HIF-1 signaling pathway',
            '04068':'FoxO signaling pathway',
            '04020':'Calcium signaling pathway',
            '04070':'Phosphatidylinositol signaling system',
            '04072':'Phospholipase D signaling pathway',
            '04071':'Sphingolipid signaling pathway',
            '04024':'cAMP signaling pathway',
            '04022':'cGMP-PKG signaling pathway',
            '04151':'PI3K-Akt signaling pathway',
            '04152':'AMPK signaling pathway',
            '04150':'mTOR signaling pathway',
            '04075':'Plant hormone signal transduction',
            # 3.3 Signaling molecules and interaction
            '04080':'Neuroactive ligand-receptor interaction',
            '04060':'Cytokine-cytokine receptor interaction',
            '04061':'Viral protein interaction with cytokine and cytokine receptor New!',
            '04512':'ECM-receptor interaction',
            '04514':'Cell adhesion molecules (CAMs)',
            # 4. Cellular Processes
            # 4.1 Transport and catabolism
            '04144':'Endocytosis',
            '04145':'Phagosome',
            '04142':'Lysosome',
            '04146':'Peroxisome',
            '04140':'Autophagy - animal',
            '04138':'Autophagy - yeast',
            '04136':'Autophagy - other',
            '04137':'Mitophagy - animal',
            '04139':'Mitophagy - yeast',
            # 4.2 Cell growth and death
            '04110':'Cell cycle',
            '04111':'Cell cycle - yeast',
            '04112':'Cell cycle - Caulobacter',
            '04113':'Meiosis - yeast',
            '04114':'Oocyte meiosis',
            '04210':'Apoptosis',
            '04214':'Apoptosis - fly',
            '04215':'Apoptosis - multiple species',
            '04216':'Ferroptosis',
            '04217':'Necroptosis',
            '04115':'p53 signaling pathway',
            '04218':'Cellular senescence',
            # 4.3 Cellular community - eukaryotes
            '04510':'Focal adhesion',
            '04520':'Adherens junction',
            '04530':'Tight junction',
            '04540':'Gap junction',
            '04550':'Signaling pathways regulating pluripotency of stem cells',
            # 4.4 Cellular community - prokaryotes
            '02024':'Quorum sensing',
            '05111':'Biofilm formation - Vibrio cholerae',
            '02025':'Biofilm formation - Pseudomonas aeruginosa',
            '02026':'Biofilm formation - Escherichia coli',
            # 4.5 Cell motility
            '02030':'Bacterial chemotaxis',
            '02040':'Flagellar assembly',
            '04810':'Regulation of actin cytoskeleton',
            # 5. Organismal Systems
            # 5.1 Immune system
            '04640':'Hematopoietic cell lineage',
            '04610':'Complement and coagulation cascades',
            '04611':'Platelet activation',
            '04620':'Toll-like receptor signaling pathway',
            '04624':'Toll and Imd signaling pathway',
            '04621':'NOD-like receptor signaling pathway',
            '04622':'RIG-I-like receptor signaling pathway',
            '04623':'Cytosolic DNA-sensing pathway',
            '04625':'C-type lectin receptor signaling pathway',
            '04650':'Natural killer cell mediated cytotoxicity',
            '04612':'Antigen processing and presentation',
            '04660':'T cell receptor signaling pathway',
            '04658':'Th1 and Th2 cell differentiation',
            '04659':'Th17 cell differentiation',
            '04657':'IL-17 signaling pathway',
            '04662':'B cell receptor signaling pathway',
            '04664':'Fc epsilon RI signaling pathway',
            '04666':'Fc gamma R-mediated phagocytosis',
            '04670':'Leukocyte transendothelial migration',
            '04672':'Intestinal immune network for IgA production',
            '04062':'Chemokine signaling pathway',
            # 5.2 Endocrine system
            '04911':'Insulin secretion',
            '04910':'Insulin signaling pathway',
            '04922':'Glucagon signaling pathway',
            '04923':'Regulation of lipolysis in adipocytes',
            '04920':'Adipocytokine signaling pathway',
            '03320':'PPAR signaling pathway',
            '04912':'GnRH signaling pathway',
            '04913':'Ovarian steroidogenesis',
            '04915':'Estrogen signaling pathway',
            '04914':'Progesterone-mediated oocyte maturation',
            '04917':'Prolactin signaling pathway',
            '04921':'Oxytocin signaling pathway',
            '04926':'Relaxin signaling pathway',
            '04918':'Thyroid hormone synthesis',
            '04919':'Thyroid hormone signaling pathway',
            '04928':'Parathyroid hormone synthesis, secretion and action',
            '04916':'Melanogenesis',
            '04924':'Renin secretion',
            '04614':'Renin-angiotensin system',
            '04925':'Aldosterone synthesis and secretion',
            '04927':'Cortisol synthesis and secretion',
            # 5.3 Circulatory system
            '04260':'Cardiac muscle contraction',
            '04261':'Adrenergic signaling in cardiomyocytes',
            '04270':'Vascular smooth muscle contraction',
            # 5.4 Digestive system
            '04970':'Salivary secretion',
            '04971':'Gastric acid secretion',
            '04972':'Pancreatic secretion',
            '04976':'Bile secretion',
            '04973':'Carbohydrate digestion and absorption',
            '04974':'Protein digestion and absorption',
            '04975':'Fat digestion and absorption',
            '04979':'Cholesterol metabolism',
            '04977':'Vitamin digestion and absorption',
            '04978':'Mineral absorption',
            # 5.5 Excretory system
            '04962':'Vasopressin-regulated water reabsorption',
            '04960':'Aldosterone-regulated sodium reabsorption',
            '04961':'Endocrine and other factor-regulated calcium reabsorption',
            '04964':'Proximal tubule bicarbonate reclamation',
            '04966':'Collecting duct acid secretion',
            # 5.6 Nervous system
            '04724':'Glutamatergic synapse',
            '04727':'GABAergic synapse',
            '04725':'Cholinergic synapse',
            '04728':'Dopaminergic synapse',
            '04726':'Serotonergic synapse',
            '04720':'Long-term potentiation',
            '04730':'Long-term depression',
            '04723':'Retrograde endocannabinoid signaling',
            '04721':'Synaptic vesicle cycle',
            '04722':'Neurotrophin signaling pathway',
            # 5.7 Sensory system
            '04744':'Phototransduction',
            '04745':'Phototransduction - fly',
            '04740':'Olfactory transduction',
            '04742':'Taste transduction',
            '04750':'Inflammatory mediator regulation of TRP channels',
            # 5.8 Development and regeneration
            '04320':'Dorso-ventral axis formation',
            '04360':'Axon guidance',
            '04361':'Axon regeneration New!',
            '04380':'Osteoclast differentiation',
            # 5.9 Aging
            '04211':'Longevity regulating pathway',
            '04212':'Longevity regulating pathway - worm',
            '04213':'Longevity regulating pathway - multiple species',
            # 5.10 Environmental adaptation
            '04710':'Circadian rhythm',
            '04713':'Circadian entrainment',
            '04711':'Circadian rhythm - fly',
            '04712':'Circadian rhythm - plant',
            '04714':'Thermogenesis',
            '04626':'Plant-pathogen interaction',
            # 6. Human Diseases
            # 6.1 Cancer: overview
            '05200':'Pathways in cancer',
            '05202':'Transcriptional misregulation in cancer',
            '05206':'MicroRNAs in cancer',
            '05205':'Proteoglycans in cancer',
            '05204':'Chemical carcinogenesis',
            '05203':'Viral carcinogenesis',
            '05230':'Central carbon metabolism in cancer',
            '05231':'Choline metabolism in cancer',
            '05235':'PD-L1 expression and PD-1 checkpoint pathway in cancer',
            # 6.2 Cancer: specific types
            '05210':'Colorectal cancer',
            '05212':'Pancreatic cancer',
            '05225':'Hepatocellular carcinoma',
            '05226':'Gastric cancer',
            '05214':'Glioma',
            '05216':'Thyroid cancer',
            '05221':'Acute myeloid leukemia',
            '05220':'Chronic myeloid leukemia',
            '05217':'Basal cell carcinoma',
            '05218':'Melanoma',
            '05211':'Renal cell carcinoma',
            '05219':'Bladder cancer',
            '05215':'Prostate cancer',
            '05213':'Endometrial cancer',
            '05224':'Breast cancer',
            '05222':'Small cell lung cancer',
            '05223':'Non-small cell lung cancer',
            # 6.3 Immune disease
            '05310':'Asthma',
            '05322':'Systemic lupus erythematosus',
            '05323':'Rheumatoid arthritis',
            '05320':'Autoimmune thyroid disease',
            '05321':'Inflammatory bowel disease (IBD)',
            '05330':'Allograft rejection',
            '05332':'Graft-versus-host disease',
            '05340':'Primary immunodeficiency',
            # 6.4 Neurodegenerative disease
            '05010':'Alzheimer disease',
            '05012':'Parkinson disease',
            '05014':'Amyotrophic lateral sclerosis (ALS)',
            '05016':'Huntington disease',
            '05020':'Prion diseases',
            # 6.5 Substance dependence
            '05030':'Cocaine addiction',
            '05031':'Amphetamine addiction',
            '05032':'Morphine addiction',
            '05033':'Nicotine addiction',
            '05034':'Alcoholism',
            # 6.6 Cardiovascular disease
            '05418':'Fluid shear stress and atherosclerosis',
            '05410':'Hypertrophic cardiomyopathy (HCM)',
            '05412':'Arrhythmogenic right ventricular cardiomyopathy (ARVC)',
            '05414':'Dilated cardiomyopathy (DCM)',
            '05416':'Viral myocarditis',
            # 6.7 Endocrine and metabolic disease
            '04930':'Type II diabetes mellitus',
            '04940':'Type I diabetes mellitus',
            '04950':'Maturity onset diabetes of the young',
            '04932':'Non-alcoholic fatty liver disease (NAFLD)',
            '04931':'Insulin resistance',
            '04933':'AGE-RAGE signaling pathway in diabetic complications',
            '04934':'Cushing syndrome',
            # 6.8 Infectious disease: bacterial
            '05110':'Vibrio cholerae infection',
            '05120':'Epithelial cell signaling in Helicobacter pylori infection',
            '05130':'Pathogenic Escherichia coli infection',
            '05132':'Salmonella infection',
            '05131':'Shigellosis',
            '05135':'Yersinia infection New!',
            '05133':'Pertussis',
            '05134':'Legionellosis',
            '05150':'Staphylococcus aureus infection',
            '05152':'Tuberculosis',
            '05100':'Bacterial invasion of epithelial cells',
            # 6.9 Infectious disease: viral
            '05166':'Human T-cell leukemia virus 1 infection',
            '05170':'Human immunodeficiency virus 1 infection',
            '05162':'Measles',
            '05164':'Influenza A',
            '05161':'Hepatitis B',
            '05160':'Hepatitis C',
            '05168':'Herpes simplex virus 1 infection',
            '05163':'Human cytomegalovirus infection',
            '05167':'Kaposi sarcoma-associated herpesvirus infection',
            '05169':'Epstein-Barr virus infection',
            '05165':'Human papillomavirus infection',
            # 6.10 Infectious disease: parasitic
            '05146':'Amoebiasis',
            '05144':'Malaria',
            '05145':'Toxoplasmosis',
            '05140':'Leishmaniasis',
            '05142':'Chagas disease (American trypanosomiasis)',
            '05143':'African trypanosomiasis',
            # 6.11 Drug resistance: antimicrobial
            '01501':'beta-Lactam resistance',
            '01502':'Vancomycin resistance',
            '01503':'Cationic antimicrobial peptide (CAMP) resistance',
            # 6.12 Drug resistance: antineoplastic
            '01521':'EGFR tyrosine kinase inhibitor resistance',
            '01524':'Platinum drug resistance',
            '01523':'Antifolate resistance',
            '01522':'Endocrine resistance',
            # 7. Drug Development
            # 7.1 Chronology: Antiinfectives
            '07011':'Penicillins',
            '07012':'Cephalosporins - parenteral agents',
            '07013':'Cephalosporins - oral agents',
            '07021':'Aminoglycosides',
            '07019':'Tetracyclines',
            '07020':'Macrolides and ketolides',
            '07014':'Quinolones',
            '07023':'Rifamycins',
            '07026':'Antifungal agents',
            '07044':'Antiviral agents',
            '07053':'Anti-HIV agents',
            # 7.2 Chronology: Antineoplastics
            '07040':'Antineoplastics - alkylating agents',
            '07041':'Antineoplastics - antimetabolic agents',
            '07042':'Antineoplastics - agents from natural products',
            '07043':'Antineoplastics - hormones',
            '07045':'Antineoplastics - protein kinases inhibitors',
            # 7.3 Chronology: Nervous system agents
            '07032':'Hypnotics',
            '07030':'Anxiolytics',
            '07033':'Anticonvulsants',
            '07015':'Local analgesics',
            '07039':'Opioid analgesics',
            '07028':'Antipsychotics',
            '07029':'Antipsychotics - phenothiazines',
            '07031':'Antipsychotics - butyrophenones',
            '07027':'Antidepressants',
            '07056':'Agents for Alzheimer-type dementia',
            '07057':'Antiparkinsonian agents',
            # 7.4 Chronology: Other drugs
            '07055':'Sulfonamide derivatives - overview',
            '07016':'Sulfonamide derivatives - sulfa drugs',
            '07017':'Sulfonamide derivatives - diuretics',
            '07018':'Sulfonamide derivatives - hypoglycemic agents',
            '07037':'Antiarrhythmic drugs',
            '07038':'Antiulcer drugs',
            '07046':'Immunosuppressive agents',
            '07047':'Osteoporosis drugs',
            '07048':'Antimigraines',
            '07049':'Antithrombosis agents',
            '07050':'Antirheumatics - DMARDs and biological agents',
            '07051':'Antidiabetics',
            '07052':'Antidyslipidemic agents',
            '07054':'Antiglaucoma agents',
            # 7.5 Target-based classification: G protein-coupled receptors
            '07220':'Cholinergic and anticholinergic drugs',
            '07215':'alpha-Adrenergic receptor agonists/antagonists',
            '07214':'beta-Adrenergic receptor agonists/antagonists',
            '07213':'Dopamine receptor agonists/antagonists',
            '07212':'Histamine H1 receptor antagonists',
            '07227':'Histamine H2/H3 receptor agonists/antagonists',
            '07211':'Serotonin receptor agonists/antagonists',
            '07228':'Eicosanoid receptor agonists/antagonists',
            '07224':'Opioid receptor agonists/antagonists',
            '07229':'Angiotensin receptor and endothelin receptor antagonists',
            # 7.6 Target-based classification: Nuclear receptors
            '07225':'Glucocorticoid and mineralocorticoid receptor agonists/antagonists',
            '07226':'Progesterone, androgen and estrogen receptor agonists/antagonists',
            '07223':'Retinoic acid receptor (RAR) and retinoid X receptor (RXR) agonists/antagonists',
            '07222':'Peroxisome proliferator-activated receptor (PPAR) agonists',
            # 7.7 Target-based classification: Ion channels
            '07221':'Nicotinic cholinergic receptor antagonists',
            '07230':'GABA-A receptor agonists/antagonists',
            '07036':'Calcium channel blocking drugs',
            '07231':'Sodium channel blocking drugs',
            '07232':'Potassium channel blocking and opening drugs',
            '07235':'N-Metyl-D-aspartic acid receptor antagonists',
            # 7.8 Target-based classification: Transporters
            '07233':'Ion transporter inhibitors',
            '07234':'Neurotransmitter transporter inhibitors',
            # 7.9 Target-based classification: Enzymes
            '07216':'Catecholamine transferase inhibitors',
            '07219':'Cyclooxygenase inhibitors',
            '07024':'HMG-CoA reductase inhibitors',
            '07217':'Renin-angiotensin system inhibitors',
            '07218':'HIV protease inhibitors',
            # 7.10 Structure-based classification
            '07025':'Quinolines',
            '07034':'Eicosanoids',
            '07035':'Prostaglandins',
            # 7.11 Skeleton-based classification
            '07110':'Benzoic acid family',
            '07112':'1,2-Diphenyl substitution family',
            '07114':'Naphthalene family',
            '07117':'Benzodiazepine family'}
    
        self.default_maps = ['00010','00020','00030','00040','00051','00052','00053',
            '00500','00520','00620','00630','00640','00650','00660','00562','00190',
            '00195','00196','00710','00720','00680','00910','00920','00061','00062',
            '00071','00072','00073','00100','00120','00121','00140','00561','00564',
            '00565','00600','00590','00591','00592','01040','00230','00240','00250',
            '00260','00270','00280','00290','00310','00220','00330','00340','00350',
            '00360','00380','00400','00410','00430','00440','00450','00460','00471',
            '00472','00473','00480','00510','00513','00512','00515','00514','00532',
            '00534','00533','00531','00563','00601','00603','00604','00540','00550',
            '00511','00571','00572','00730','00740','00750','00760','00770','00780',
            '00785','00790','00670','00830','00860','00130','00900','00902','00909',
            '00904','00906','00905','00981','00908','00903','00281','01052','00522',
            '01051','01059','01056','01057','00253','00523','01054','01053','01055',
            '00940','00945','00941','00944','00942','00943','00901','00403','00950',
            '00960','01058','00232','00965','00966','00402','00311','00332','00261',
            '00331','00521','00524','00525','00231','00401','00404','00405','00333',
            '00254','00998','00999','00362','00627','00364','00625','00361','00623',
            '00622','00633','00642','00643','00791','00930','00363','00621','00626',
            '00624','00365','00984','00980','00098','00098','01010','01060','01061',
            '01062','01063','01064','01065','01066','01070','03020','03022','03040',
            '03010','00970','03013','03015','03008','03060','04141','04130','04120',
            '04122','03050','03018','03030','03410','03420','03430','03440','03450',
            '03460','02010','02060','03070','02020','04014','04015','04010','04011',
            '04012','04310','04330','04340','04350','04390','04392','04370','04371',
            '04630','04064','04668','04066','04068','04020','04070','04072','04071',
            '04024','04022','04151','04152','04150','04075','04080','04060','04061',
            '04512','04514','04144','04145','04142','04146','04138','04136','04139',
            '04110','04111','04112','04113','04114','04210','04215','04216','04217',
            '04115','04218','04510','04520','04530','04540','04550','02024','05111',
            '02025','02026','02030','02040','04810']
        
    def keggid2ko(self, kegg_ids, step = 150):
        '''
        Converts KEGG_ID genes to Ortholog KO ID from KEGG
        :param KEGG_ID: (list) - KEGG ID genes
        :param step: (int) - will convert "step" KEGG IDs at a time
        :return: (list) - (list,list) - KEGG ID genes converted and ko IDs
        '''
        result = list()
        i = 0
        print('Converting {:d} KEGG IDs to KOs through the KEGG API.'.format(len(kegg_ids)))
        pbar = ProgressBar()
        for i in pbar(range(0, len(kegg_ids) - step, step)):
            try:
                result += kegg_link("ko", kegg_ids[i:i+step]).read().split("\n")[:-1]
            except:
                print('Broke at index ' + str(i))
                result = [relation.split('\t') for relation in result]
                return list(map(list, zip(*result)))
            
        result += kegg_link("ko", kegg_ids[len(kegg_ids) - step:]).read().split("\n")[:-1]; open('conversion.txt', 'w').write('\n'.join(result))
        result = [relation.split('\t') for relation in result]
        return list(map(list, zip(*result)))

    def keggid2ko_mapping(self, kegg_ids, step = 150, output_file = None):
        '''
        Adaptation from keggid2ko to include writing of mapping if it has already
        been performed
        :param KEGG_ID: (list) - KEGG ID genes
        :param step: (int) - will convert "step" KEGG IDs at a time
        :param output_file: (str) - Name of file to output mapping to
        :return: (list) - (list,list) - KEGG ID genes converted and ko IDs
        '''
        if output_file and os.path.isfile(output_file):
            return pd.read_csv(output_file, sep = '\t')
        mapping = pd.DataFrame(self.keggid2ko(kegg_ids, step = step)).transpose()
        mapping.columns = ['KEGG ID', 'KO']
        if output_file:
            mapping.to_csv(output_file, sep = '\t', index = False)
        return mapping
    
    def ko2ec(self, ko_ID):
        '''
        Converts kegg ortholog id to EC numers
        :param ko_ID: list of kegg ortholog ids
        :return: dic associating ortholog kegg id with list
        of assotiated EC numbers
        '''
        ko_dic = {}
        conversion = kegg_link("enzyme", ko_ID).read().split("\n")[:-1]
        ko_dic = {result.split("\t")[0].strip("ko:"):result.split("\t")[1].upper()
                for result in conversion if len(result) > 1}                    # if len(result) > 1 -> some KOs might not have EC numbers
        ko_dic = {}
        for result in conversion:
            if len(result) > 1:
                ko_dic[result.split("\t")[0].strip("ko:")] = result.split("\t")[1].upper()
        return ko_dic


kp = KeggPathway()

class KeggMap():
    '''
    This class 
    '''

    def __init__(self, pathway_ID, **kwargs):
        '''
        Initialize object
        :param pathway_ID: (str) - KEGG Pathway ID
        '''
        self.__dict__ = kwargs
        self.pathway_ID = pathway_ID[-5:] if len(pathway_ID) > 5 else pathway_ID
        self.set_pathway(pathway_ID)

    ############################################################################
    ####                              Helper                                ####
    ############################################################################
    
    def set_bgcolor(self, pathway_element, color):
        '''
        Sets graphic element background color
        :param pathway_element: kegg pathway xml element object
        :param color: color to be used in rgb
        '''
        pathway_element.graphics[0].bgcolor = color

    def set_fgcolor(self, pathway_element, color):
        '''
        Sets graphic element contour color
        :param pathway_element: kegg pathway xml element object
        :param color:  color to be used in rgb
        '''
        pathway_element.graphics[0].fgcolor = color
        
    def set_colors(self, colors = [], ncolor = 1):
        '''
        Creates list of hex colors to be used, using matplotlib or using custom colors
        :param colors: list of hex colors
        :param ncolor: int indicating the amount of hex colors that should be created
        :return: returns list with hex color codes
        '''
        if len(colors) == 0:
            # if no colors are given creates a list of hex colors with ncolor
            # ammount from matplotlib discrete colormaps. If ncolor > 12
            # a continuous colormap is used instead.
            if ncolor <= 8:
                pastel2 = cm.get_cmap('Pastel2', 8)
                return [to_hex(pastel2(i)) for i in range(ncolor)]
            elif ncolor <= 12:
                set3 = cm.get_cmap("Set3", 12)
                return [to_hex(set3(i)) for i in range(ncolor)]
            else:
                rainbow = cm.get_cmap("rainbow",ncolor)
                return [to_hex(rainbow(i)) for i in range(ncolor)]
        else:
            # validates hex values and returns the original list
            isvalidhex = True
            i = 0
            while isvalidhex:
                match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', str(colors[i]))
                if not match:
                    isvalidhex = False
                i += 1
            if isvalidhex:
                return colors
            else:
                raise Exception("Colors aren't valid hex codes")

    def conv_value_rgb(self, value, colormap, norm):
        '''
        Normalizes values in a vector and transforms it into corresponding
        hex color using given colormap
        :param value: numpy vector from dataframe apply function with expression
        values
        :param colormap: matplotlib colormap object
        :param norm: matplotlib Normalize or LogNorm object
        :return: returns hex colors in vector format
        '''
        return value.apply(norm).apply(colormap)

    def conv_rgb_hex(self, rgb):
        '''
        converts rgb into hex color code in vector
        :param rgb: rgb value in vector resulting from pandas dataframe apply
        :return: vector with converted hex color code
        '''
        return rgb.apply(to_hex)

    def load_pathway(self, pathway_ID, organism_ID = ""):
        '''
        Downloads pathway kgml from KEGG and reads it
        :param pathway_ID: (str) - suffix Pathway ID, int part
        :param organism_ID: (str) - preffix Pathway ID, str part
        :return: (object) - pathway kgml parsed
        '''
        if not organism_ID:
            pathway = kegg_get("ko" + pathway_ID, "kgml")
            return KGML_parser.read(pathway)
        else:
            try:
                pathway = kegg_get(organism_ID + pathway_ID, "kgml")
                return KGML_parser.read(pathway)
            except:
                print("Invalid IDs")
                
    def reset_pathway(self):
        '''
        Resets pathway state
        '''
        self.set_pathway(self.pathway_ID)

    def organismo_genes(self, kegg_ID):
        '''
        returns list of organism ids from a given list of kegg gene id
        :param kegg_ID:list with kegg gene id
        :return:list with kegg org ids
        '''
        organisms = []
        for gene_ID in kegg_ID:
            organisms.append(gene_ID.split(":")[0])
        return list(set(organisms))

    def ortholog_dic(self):
        '''
        Dictionary associating all ortholog names in pathway element
        with its first ortholog
        :return: dic associaiting orthologs present in pathway elements
        '''
        orthologs_dic = {}
        for ortholog_rec in self.pathway.orthologs:
            ortholog_name = ortholog_rec.name.split(" ")[0]
            for ortholog in ortholog_rec.name.split(" "):
                if ortholog not in orthologs_dic.keys():
                    orthologs_dic[ortholog] = ortholog_name
        return orthologs_dic
    
    ############################################################################
    ####                            Sets                                    ####
    ############################################################################

    def set_pathway(self, pathway_ID):
        '''
        Set pathway with Kegg Pathway ID
        :param pathway_ID: (str) Kegg Pathway ID
        '''
        self.pathway = self.load_pathway(self.pathway_ID)                       # get the KGML
        ko = []
        self.ko_boxes = {}
        for i in range(len(self.pathway.orthologs)):
            self.set_bgcolor(self.pathway.orthologs[i], "#ffffff")              # set all boxes to white
            self.set_fgcolor(self.pathway.orthologs[i], "#ffffff")              # ditto
            orthologs_in_box = [ide[3:] for ide in self.pathway.orthologs[i].name.split()]  # 'ko:K16157 ko:K16158 ko:K16159' -> ['K16157', 'K16158', 'K16159']
            for ortholog in orthologs_in_box:
                if ortholog not in self.ko_boxes.keys():
                    self.ko_boxes[ortholog] = list()
                self.ko_boxes[ortholog].append(i)                               # {'K16157':[0,13,432], 'K16158':[4,13,545]}
            ko.append(self.pathway.orthologs[i].graphics[0].name.strip("."))    # 'K16157...' -> 'K16157'
        
        # Set text in boxes to EC numbers
        ko_to_ec = kp.ko2ec(ko)                                                 # {'K16157':'ec:1.14.13.25'}
        for ortholog_rec in self.pathway.orthologs:
            ko = ortholog_rec.graphics[0].name.strip(".")
            if ko in ko_to_ec.keys():
                ortholog_rec.graphics[0].name = ko_to_ec[ko]

    ############################################################################
    ####                    Graphical Manipulation                          ####
    ############################################################################

    def create_tile_box(self, record):
        '''
        Create box graphical element in pathway to draw the box countour and
        give the correct name
        :param record: graphical element to be created
        '''
        newrecord = KGML_pathway.Graphics(record)
        newrecord.name = record.graphics[0].name
        newrecord.type = "rectangle"
        newrecord.width = record.graphics[0].width
        newrecord.height = record.graphics[0].height
        newrecord.y = record.graphics[0].y
        newrecord.x = record.graphics[0].x
        newrecord.bgcolor = "#FFFFFF00"
        newrecord.fgcolor = "#000000"
        record.graphics.append(newrecord)
        record.graphics[0].bgcolor = "#FFFFFF00"
        record.graphics[0].fgcolor = "#FFFFFF00"
        record.graphics[0].name = ""

    def create_box_heatmap(self, rec_old, nrboxes, i, paired = True):
        '''
        Helper function for creating heatmap, draws one expression value in its
        correct position on the bigger parent box
        :param rec_old: graphical element object to be used as reference
        :param nrboxes: int nr of boxes to be drawed
        :param i: int internal number of movements of the box given by the for loop
        :return: graphical element object
        '''
        movement_steps = rec_old.graphics[0].width / (nrboxes * (2 if paired else 1))
        newrecord = KGML_pathway.Graphics(rec_old)
        newrecord.name = ""
        newrecord.type = "rectangle"
        newrecord.width = movement_steps * 1.3 * (2 if paired else 1)           # sub-boxes width, adjusted by a factor that experimentally fixed well in the representations
        newrecord.height = rec_old.graphics[0].height
        newrecord.y = rec_old.graphics[0].y
        newrecord.x = (i * movement_steps) + rec_old.graphics[0].x
        newrecord.fgcolor = "#FFFFFF00"
        return newrecord

    ############################################################################
    ####                          Operations                                ####
    ############################################################################

    def pathway_pdf(self, filename, imagemap = True, orthologs = True, 
                    compounds = True, maps = True, reactions = True):
        '''
        Prints current pathway to PDF file
        :param filename: (str) - PDF filename
        :param imagemap: (bol) - Print imagemap
        :param orthologs: (bol) - Print orthologs
        :param compounds: (bol) - Print compounds
        :param maps: (bol) - Print maps
        :param reactions: (bol) - Print reactions ???
        :return: creates PDF file with current pathway
        '''
        #TODO Verificar o parametro reactions
        pathway_pdf = KGMLCanvas(self.pathway)
        pathway_pdf.import_imagemap = imagemap
        pathway_pdf.label_orthologs = orthologs
        pathway_pdf.label_compounds = compounds
        pathway_pdf.label_maps = maps
        pathway_pdf.label_reaction_entries = reactions
        pathway_pdf.draw(filename)

    def taxon_map(self, organism_ID, color = ""):
        '''
        Paints a map with the boxes for which the given taxon has KOs identified
        :param organism_ID: KEGG ID of the taxon
        :param color: hex color code
        '''
        #TODO adaptar para recolher varios organismos

        color1 = self.set_colors(color)[0]

        pathway_organism = self.load_pathway(self.pathway_ID, organism_ID)
        organism_kegg_ID = []
        for gene_rec in pathway_organism.genes:
            for gene in gene_rec.name.split(" "):
                organism_kegg_ID.append(gene)
        organism_ortholog_ID = kp.keggid2ko_mapping(
                organism_kegg_ID, output_file = self.output + '/kegg2ko.tsv')

        for ortholog_rec in self.pathway.orthologs:
            for ortholog in ortholog_rec.name.split(" "):
                if ortholog in organism_ortholog_ID['KO']:
                    self.set_bgcolor(ortholog_rec, color1)
                    self.set_fgcolor(ortholog_rec, color1)

    def pathway_genes_organismo(self, KEGG_ID, maxshared = 5, color = []):
        '''
        Given a list of gene_kegg ids colors the the pathway map with
        a specific color for organism
        :param KEGG_ID: list of kegg gene id
        :param maxshared: naximum number of boxes to be drawed
        :param color: list of hex color codes to be used to color the pathway map
        '''
        coverted_KEGG_ID, orthologs_ID = kp.keggid2ko_mapping(
                KEGG_ID, output_file = self.output + '/kegg2ko.tsv')
        organisms = self.organismo_genes(coverted_KEGG_ID)

        # Set colors
        organism_colors = {}
        if len(color) < len(organisms):
            colors = self.set_colors([], len(organisms))
        else:
            colors = self.set_colors(color)
        i = 0
        for organism in organisms:
            if organism not in organism_colors.keys():
                organism_colors[organism] = colors[i]
                i += 1

        # Dictionario Ortholog - [organismo]
        ortholog_dic = {}
        for i in range(len(orthologs_ID)):
            if orthologs_ID[i-1] in ortholog_dic.keys():
                ortholog_dic[orthologs_ID[i - 1]].append(organisms[i - 1])
            else:
                ortholog_dic[orthologs_ID[i - 1]] = [organisms[i - 1]]

        for ortholog_rec in self.pathway.orthologs:
            # organismos por ortholog_rec
            organisms = []
            for ortholog in ortholog_rec.name.split(" "):
                if ortholog in ortholog_dic.keys():
                    organisms += ortholog_dic[ortholog]
            # Se tiver organismos para o ortholog_rec
            if len(organisms) > 0:
                if len(organisms) <= maxshared:
                    organism_size = len(organisms)
                else:
                    organism_size = maxshared
                organisms.sort()
                if organism_size%2==0:
                    n = 0
                    for i in range(-organism_size+1, organism_size,2):
                        newrecord = self.create_pair_boxes(ortholog_rec, organism_size, i)
                        newrecord.bgcolor = organism_colors[organisms[n]]
                        ortholog_rec.graphics.append(newrecord)
                        n+=1
                    self.create_tile_box(ortholog_rec)

                if organism_size% 2 != 0:
                    n= 0
                    for i in range(int(-(organism_size-1)/2),int((organism_size+1)/2)):
                        newrecord = self.create_odd_boxes(ortholog_rec, organism_size, i)
                        newrecord.bgcolor = organism_colors[organisms[n]]
                        ortholog_rec.graphics.append(newrecord)
                        n += 1
                    self.create_tile_box(ortholog_rec)

    def pathway_box_list(self, dic_box_items, items, maxshared = 10, color = []):
        '''
        Represents items in the pathway map
        :param dic_box_items: dic with keys correspoding to int index of the list
        of ortholog elements in pathway and values correspoding to elements
        of items
        :param items: list of items to be represented and given a specific color
        :param maxshared: int representing the maximum number of boxes to be
        represented in each element in the pathway map
        :param color: list of costum colors to be used to color the elements
        '''
        colors = self.set_colors(color, len(items))
        dic_colors = {items[i] : colors[i] for i in range(len(items))}

        for box in dic_box_items.keys():
            number_types = len(dic_box_items[box])
            if number_types > maxshared:
                number_types = maxshared
            if number_types % 2 == 0:
                n = 0
                for i in range(-number_types + 1, number_types, 2):
                    newrecord = self.create_pair_boxes(self.pathway.orthologs[box], number_types, i)
                    newrecord.bgcolor = dic_colors[dic_box_items[box][n]]
                    self.pathway.orthologs[box].graphics.append(newrecord)
                    n += 1
                self.create_tile_box(self.pathway.orthologs[box])

            else:
                n = 0
                for i in range(int(-(number_types - 1) / 2),
                               int((number_types + 1) / 2)):
                    newrecord = self.create_odd_boxes(self.pathway.orthologs[box], number_types, i)
                    newrecord.bgcolor = dic_colors[dic_box_items[box][n]]
                    self.pathway.orthologs[box].graphics.append(newrecord)
                    n += 1
                self.create_tile_box(self.pathway.orthologs[box])

    def pathway_boxes_diferential(self, dataframe, log = False, colormap = "coolwarm"):
        '''
        Represents expression values present in a dataframe in the
        pathway map
        :param dataframe: pandas DataFrame with each column representing a sample
        and index corresponding to int list index of the ortholog elment in the
        pathway
        :param log: bol providing the option for a log normalization of data
        :param colormap: str representing a costum matplotlib colormap to be used
        '''
        
        if log:
            norm = cm.colors.LogNorm(vmin=dataframe.min().min(), vmax=dataframe.max().max())
        else:
            norm = cm.colors.Normalize(vmin=dataframe.min().min(), vmax=dataframe.max().max())

        colormap = cm.get_cmap(colormap)
        dataframe = dataframe.apply(self.conv_value_rgb, args=(colormap, norm))
        dataframe = dataframe.apply(self.conv_rgb_hex)
        
        dataframe = dataframe[dataframe.columns.tolist()]
        
        nrboxes = len(dataframe.columns.tolist())                               # number of samples

        for box in dataframe.index.tolist():
            colors = dataframe.loc[box].tolist()
            paired = True if nrboxes % 2 == 0 else False
            for i in range(nrboxes):
                newrecord = self.create_box_heatmap(self.pathway.orthologs[box], nrboxes, 
                                                    i * 2 - (nrboxes - 1) if paired     # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7
                                                    else i - int(nrboxes / 2),          # if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                                                    paired = paired)
                newrecord.bgcolor = colors[i]
                self.pathway.orthologs[box].graphics.append(newrecord)
            self.create_tile_box(self.pathway.orthologs[box])

    def grey_boxes(self, box_list):
        for i in box_list:
            self.set_bgcolor(self.pathway.orthologs[i], "#7c7272")
            self.set_fgcolor(self.pathway.orthologs[i], "#7c7272")


################################################################################
################################################################################

class MoscaData():
    def __init__(self, file = None, pathway=[], **kwargs):
        self.__dict__ = kwargs
        if file:
            self.data = self.load_file(file)                                    #            Protein ID  ...               KEGG ID
                                                                                # KO                     ...
                                                                                # ko:K01126  A0A0P0GQB3  ...  bcel:BcellWH2_02985;
        if len(pathway) > 0:
            self.pathway = pathway
            
    # TODO - whats happens when no pathway is specified
    # handle MOSCA output/load_file

    def load_file(self, file):
        '''
        Load MOSCA .csv file and prepares it for operations
        :param file: .csv file from MOSCA analysis
        return: pandas DataFrame object containing the formated .csv data
        '''
        df = pd.read_csv(file, sep = '\t')
        df = df[df['Cross-reference (KEGG)'].notnull()]
        data = kp.keggid2ko_mapping(df['Cross-reference (KEGG)'].tolist(),      #          KEGG ID         KO
                                    output_file = self.output + '/kegg2ko.tsv') # 0  sfu:Sfum_1080  ko:K03738
        data['KEGG ID'] += ';'
        df = pd.merge(df, data, left_on = 'Cross-reference (KEGG)', right_on = 'KEGG ID')   # inner join
        df.index = [ide.split(':')[1] for ide in df['KO']]                      # ko:K00001 -> K00001
        return df

    def all_kegg_maps(self):
        '''
        Creates a dic with all specific kegg pathway maps and their description
        :return: pandas.DataFrame with Map ID as index and maps names as
        sole column
        '''
        maps = pd.read_csv(StringIO(kegg_list("pathway").read()), sep='\t',
                           names = ['ID', 'Name'])
        
        maps['ID'] = maps['ID'].apply(lambda x: x.split(':')[1])
        maps.set_index('ID', inplace = True)
        
        return maps

    def top_gemus(self, samples, n=10):
        '''
        Calculates top genus from samples
        :param samples: list of samples
        :param n: number of top genus to return
        :return: list of top genus
        '''
        data = self.data.groupby("Taxonomic lineage (GENUS)")[samples].sum()
        data["sums"] = data.sum(axis=1)
        data = data.sort_values(by=["sums"], ascending=False)
        if n > len(data.index.tolist()):
            n = len(data.index.tolist())
        return data.index.tolist()[:n]

    ############################################################################
    ####                              Helper                                ####
    ############################################################################

    def set_colors(self, colors = [], ncolor = 1):
        '''
        Creates list of hex colors to be used,
        using matplotlib or using costum colors
        :param colors: list of hex colors
        :param ncolor: int indicating the ammount of hex colors should be created
        :return: returns list with hex color codes
        '''
        ret_color = []
        if len(colors) == 0:
            # if no colors are given creates a list of hex colors with ncolor
            # ammount from matplotlib discrete colormaps. If ncolor > 12
            # a continuous colormap is used instead.
            if ncolor <= 8:
                pastel2 = cm.get_cmap('Pastel2', 8)
                for i in range(ncolor):
                    ret_color.append(to_hex(pastel2(i)))
                return ret_color
            elif ncolor <= 12:
                set3 = cm.get_cmap("Set3", 12)
                for i in range(ncolor):
                    ret_color.append(to_hex(set3(i)))
                return ret_color
            else:
                rainbow = cm.get_cmap("rainbow",ncolor)
                for i in range(ncolor):
                    ret_color.append(to_hex(rainbow(i)))
                return ret_color
        else:
            # validates hex values and returns the original list
            isvalidhex = True
            for hexvalue in colors:
                match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', str(hexvalue))
                if not match:
                    isvalidhex = False
            if isvalidhex:
                return colors
            else:
                raise Exception("Colors aren't valid hex codes")

    def load_pathway(self, pathway):
        '''
        Loads a given pathway map in Pathway object
        :param pathway: str of kegg map pathway
        :return: Pathway map object
        '''
        return KeggMap(pathway)

    ############################################################################
    ####                          Operations                                ####
    ############################################################################

    def genomic_potential_genus(self, samples, genera = [], nr = 10):
        '''
        Represents the genomic potential of the dataset, by coloring each
        genus with an unique color. Only the orthologs with at least one
        sample expression level above the threshold will be represented
        :param samples:list of str column names of the dataset correspoding to
        expression values
        :param genera: list of genus to represent
        :param threshold: int representing the expression threshold
        '''
        #KO for each genus type
        if len(genera) > 0:
            genera = list(set(genera))
        else:
            genera = self.top_gemus(samples, nr)
        colors = self.set_colors([], len(genera))
        for metabolic_map in self.pathway:
            try:
                pathway = self.load_pathway(metabolic_map)
                dic_boxes = {}
                present_genus = []
                for genus in genera:
                    df = self.data[self.data["Taxonomic lineage (GENUS)"] == genus][samples]
                    df = df[df.any(axis=1)]
                    orthologs = df.index.tolist()
                    for ortholog in orthologs:
                        if ortholog in pathway.ko_boxes.keys():
                            if genus not in present_genus:
                                present_genus.append(genus)
                            for box in pathway.ko_boxes[ortholog]:
                                if box in dic_boxes.keys():
                                    dic_boxes[box].append(genus)
                                else:
                                    dic_boxes[box] = [genus]
                
                pathway.pathway_box_list(dic_boxes, present_genus, len(genus), colors)
                
                df = self.data[samples]
                df = df[df.any(axis=1)]
                
                orthologs = df.index.tolist()
                grey_boxes = []
                for ortholog in orthologs:
                    if ortholog in pathway.ko_boxes.keys():
                        for box in pathway.ko_boxes[ortholog]:
                            if box not in dic_boxes.keys() and box not in grey_boxes:
                                grey_boxes.append(box)
                pathway.grey_boxes(grey_boxes)

                name_pdf = metabolic_map + "_genomic_potential"
                pathway.pathway_pdf(name_pdf)
                print("Map saved to " + name_pdf + ".pdf")
                if len(grey_boxes) > 0:
                    temp_colors = colors
                    temp_colors.append("#7c7272")
                    temp_genus = genus
                    temp_genus.append("Present in Samples")
                    self.potential_legend(temp_colors, temp_genus, str(name_pdf) + ".png")
                    
                else:
                    self.potential_legend(colors, genus, str(name_pdf) + ".png")
            except:
                print(metabolic_map + " - Failed to process")

    def genomic_potential_sample(self, samples, threshold = 0):
        '''
        Represents the genomic potential of the dataset, by coloring each sample
        with an unique color. Only orthologs iwth at least one sample expression
        level above the theshold will be represented
        :param samples: ist of str column names of the dataset correspoding to
        expression values
        :param threshold: threshold: int representing the expression threshold
        '''
        dic_box_sample = {}
        df = self.data[samples]
        for sample in samples:
            orthologs = df.loc[df[sample] > threshold].index.tolist()
            for ortholog in orthologs:
                if ortholog in self.pathway.ko_boxes.keys():
                    for box in self.pathway.ko_boxes[ortholog]:
                        if box in dic_box_sample.keys():
                            if sample not in dic_box_sample[box]:
                                dic_box_sample[box].append(sample)
                        else:
                            dic_box_sample[box] = [sample]

        self.pathway.pathway_box_list(dic_box_sample, samples)


    def differential_expression_sample(self, samples, output_folder,
                                       log=False, pathways = None):
        '''
        Represents in small heatmaps the expression levels of each sample on the
        dataset present in the given pathway map. The values can be transford to
        a log10 scale
        :param samples: list of str column names of the dataset corresponding to
        expression values
        :param output_folder: string  - name of folder to store pdfs
        :param log: bol providing the option for a log normalization of data
        '''
        
        for metabolic_map in pathways:
            pathway = self.load_pathway(metabolic_map)
            data = []
            new_index = []
            df = self.data.groupby(level = 0)[samples].sum()
            df = df[df.any(axis=1)]
            pbar = ProgressBar()
            print('Mapping {} KOs to corresponding boxes in map {}.'.format(
                    len(df.index), metabolic_map))
            for ortholog in pbar(df.index.tolist()):
                if ortholog in pathway.ko_boxes.keys():
                    for box in pathway.ko_boxes[ortholog]:
                        for ko in df.index.tolist():
                            data.append(df.loc[ko].tolist())
                            new_index.append(box)
            
            new_df = pd.DataFrame(data, columns=samples, index=new_index)       #                sample1
                                                                                # 6              151800.0
            new_df = new_df[new_df.any(axis=1)]
            new_df = new_df.groupby(level=0).sum()
            
            pathway.pathway_boxes_diferential(new_df, log)
            name_pdf = '{}/{}_differential_expression{}.pdf'.format(output_folder, 
                        metabolic_map, '_log' if log else '')
            pathway.pathway_pdf(name_pdf)
                
            print("Map saved to " + name_pdf + ".pdf")
            self.differential_colorbar(new_df, name_pdf.replace(".pdf",'png'))
                
    
    def potential_legend(self, colors, labels, filename):
        colors = colors
        f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
        handles = [f("s", colors[i]) for i in range(len(colors))]
        labels = labels
        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=True)

        def export_legend(legend, filename="legend.png", expand=[-5,-5,5,5]):
            fig  = legend.figure
            fig.canvas.draw()
            bbox  = legend.get_window_extent()
            bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
            bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(filename, dpi="figure", bbox_inches=bbox)

        export_legend(legend, filename)

    def differential_colorbar(self, dataframe, filename):
        FIGSIZE = (2,3)
        mpb = plt.pcolormesh(dataframe,cmap='coolwarm')
        fig,ax = plt.subplots(figsize=FIGSIZE)
        plt.colorbar(mpb,ax=ax)
        ax.remove()
        plt.savefig(filename,bbox_inches='tight')
                
############################################################################
####                            Test                                    ####
############################################################################
    
def test_pathway_organismo():
    test_pathway = KeggMap("map00680", output = 'debugKEGG')
    test_pathway.pathway_organismo("eco")
    test_pathway.pathway_pdf("debugKEGG/test_pathway_organismo.pdf")
    test_pathway.reset_pathway()
    test_pathway.pathway_organismo("eco", color = ["#ef0000"])
    test_pathway.pathway_pdf("debugKEGG/test_pathway_organismo_vermelho.pdf")

def test_pathway_genes_organismo():
    test_pathway = KeggMap("map00680")
    kegg_IDs = ["ppg:PputGB1_2248","ppt:PPS_3154","ppx:T1E_2051","ppun:PP4_21640","ppud:DW66_3426","sme:SMc02610", "mby:MSBRM_0152"]
    test_pathway.pathway_genes_organismo(kegg_IDs)
    test_pathway.pathway_pdf("debugKEGG/test_pathway_genes_organismo.pdf")

def test_genomic_potential_genus():
    data = MoscaData("all_info_simulated.tsv",["map00680"])
    samples_mg = ["grinder-reads0.17a", "grinder-reads0.17b", "grinder-reads0.17c"]
    data.genomic_potential_genus(samples_mg, ["Pseudomonas", "Methanosarcina","Pelolinea","Thermoplasmatales"])
    data.pathway_pdf("debugKEGG/test_genomic_potential_genus.pdf")

def test_genomic_potential_sample():
    data = MoscaData("all_info_simulated.tsv", ["map00680"])
    samples_mg = ["grinder-reads0.17a", "grinder-reads0.17b","grinder-reads0.17c"]
    data.genomic_potential_sample(samples_mg)
    data.pathway_pdf("debugKEGG/test_genomic_potential_sample.pdf")


def test_differential_expression_sample():
    data = MoscaData("marosca.tsv", output = 'debugKEGG')
    #data.set_pathway(["map00680"])

    #samples_mt = ["grinder-reads1a", "grinder-reads1b", "grinder-reads1c", "grinder-reads3a", "grinder-reads3b", "grinder-reads3c"]
    samples_mt = ["grinder-reads{}{}".format(i, j) for i in ['0.17','1','3'] for j in ['a','b','c']]
    data.differential_expression_sample(samples_mt, 'debugKEGG', log = True, pathways = ["map00680"])
    data.differential_expression_sample(samples_mt, 'debugKEGG', log = False, pathways = ["map00680"])


############################################################################
####                            Main                                    ####
############################################################################

if __name__ == '__main__':
    #kegg_link("ko", "ppg:PputGB1_2248")
    ## Test Functions
    #test_pathway_organismo()
    #test_genomic_potential_genus()
    #test_genomic_potential_sample()
    test_differential_expression_sample()

'''
pbar = ProgressBar()
for sample in samples:
    print(sample)
    pbar = ProgressBar()
    for i in pbar(marosca.index):
        if marosca.loc[i][sample] != np.nan:
            marosca.loc[i][sample] *= np.random.uniform(low=0.8, high=1.2, size = 1)[0]
'''