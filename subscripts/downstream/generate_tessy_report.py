from docxtpl import DocxTemplate
import time, pandas as pd, sys, pathlib, datetime

#Import template document + paths to files with data_1 & plots
template = DocxTemplate('/mnt/home/jevgen01/nmrl/cov_analysis/downstream/tessy_report_template.docx')
if len(sys.argv) < 2:
    sys.exit('Path to mutation report is required! (mutation_report.py can be used to produce it)')

mutation_report_path = sys.argv[1] 
output_path = pathlib.Path(mutation_report_path).parent

#Generate data_1 lineage data_1 table
data = pd.read_csv(mutation_report_path)
data.fillna("0",inplace=True)

#Combine data into report
lineage_data_table=[
    {"num":1,"name":"AY.4.2 = AY.4.2 (mutations: L452R, T478K, D614G, P681R, A222V, Y145H)","count":len(data[data['lineage'].str.contains('AY.4.2')])},
    {"num":2,"name":"B.1.1.529 = B.1.1.529: See reporting protocol, report as BA sub-lineage if possible","count":len(data[(data['lineage'] == 'B.1.1.529')])},
    {"num":3,"name":"B.1.1.7 = B.1.1.7 (mutations:del69-70,del144,N501Y,A570D,D614G,P681H,T716I,S982A,D1118H)","count":len(data[data['lineage'] =='B.1.1.7'])},
    {"num":4,"name":"B.1.1.7+E484K = B.1.1.7+E484K (mutations as B.1.1.7 and additionally E484K)","count":len(data[(data['lineage'] =='B.1.1.7') & (data['E484K'] == '1.0')])},
    {"num":5,"name":"B.1.351 = B.1.351 (defined by mutations: D80A, D215G, E484K, N501Y, A701V)","count":len(data[data['lineage'] =='B.1.351'])},
    {"num":6,"name":"B.1.427/B.1.429 = B.1.427/B.1.429 (mutations: L452R, D614G)","count":len(data[data['lineage'] =='B.1.427']) + len(data[data['lineage'] =='B.1.429'])},
    {"num":7,"name":"B.1.525 = B.1.525 (mutations:E484K, D614G, Q677H)","count":len(data[data['lineage'] =='B.1.525'])},
    {"num":8,"name":"B.1.616 = B.1.616 (mutations:D215G,D614G,142del,G669S,H66D,H655Y,N1187D,Q949R,V483A,Y144V)","count":len(data[data['lineage'] =='B.1.616'])},
    {"num":9,"name":"B.1.617 = B.1.617 lineage or any sublineage of B.1.617(common mutations:D614G,L452R,P681R)","count":len(data[data['lineage'] =='B.1.617'])},
    {"num":10,"name":"B.1.617.1 = B.1.617.1 (mutations: L452R, E484Q, D614G, P681R)","count":len(data[data['lineage'] =='B.1.617.1'])},
    {"num":11,"name":"B.1.617.2 = B.1.617.2 (mutations: L452R, T478K, D614G, P681R)","count":len(data[data['lineage'].str.contains('B.1.617.2')]) + len(data[data['lineage'].str.contains('AY.')]) - len(data[data['lineage'].str.contains('AY.4.2')])},
    {"num":12,"name":"B.1.617.3 = B.1.617.3 (mutations: L452R, E484Q, D614G, P681R)","count":len(data[data['lineage'] =='B.1.617.3'])},
    {"num":13,"name":"B.1.620 = B.1.620 (mutations: S477N, E484K, D614G, P681H)","count":len(data[data['lineage'] =='B.1.620'])},
    {"num":14,"name":"B.1.621 = B.1.621 (mutations: R346K, E484K, N501Y, D614G, P681H)","count":len(data[data['lineage'] =='B.1.621'])},
    {"num":15,"name":"BA.1 = BA.1 or B.1.1.529 with mutations del69-70, ins214EPE, S371L, G496S, T547K","count":len(data.loc[(data['lineage'].str.contains('BA.1')) | (data['lineage'] == 'B.1.1.529') & (data['A67_V70delinsVI'] == "1.0") & (data['214EPEins'] == '1.0') & (data['S371L'] == "1.0") & (data['G496S'] == '1.0') & (data['T547K'] == '1.0')])},
    {"num":16,"name":"BA.2 = BA.2 or B.1.1.529 with mutations V213G, T376A, R408S","count":len(data.loc[(data['lineage'].str.contains('BA.2')) | (data['lineage'] == 'B.1.1.529') & (data['V213G'] == "1.0") & (data['T376A'] == '1.0') & (data['R408S'] == "1.0")])},
    {"num":17,"name":"BA.3 = BA.3 or B.1.1.529 with mutations del69-70, ORF1a:A3657V, ORF3a:T22V","count":len(data.loc[(data['lineage'].str.contains('BA.3')) | (data['lineage'] == 'B.1.1.529') & (data['A67_V70delinsVI'] == "1.0") & (data['A3657V'] == '1.0') & (data['T22V'] == "1.0")])},
    {"num":18,"name":"C.37 = C.37 (mutations L452Q, F490S, D614G)","count":len(data[data['lineage'] =='C.37'])},
    {"num":19,"name":"CLUSTER_5 = Denmark cluster 5 associated with mink (defined by mutations: del 69-70, Y453F, I692V, M1229I)","count":0},
    {"num":20,"name":"E484K = detected via an SNP assay specific for E484K","count":0},
    {"num":21,"name":"N501Y = detected via an SNP assay specific for N501Y","count":0},
    {"num":22,"name":"ORF1a(del3675-3677) = Variants carrying ORF1a deletion (del 3675-3677)","count":0},
    {"num":23,"name":"P.1 = P.1 variants (L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, H655Y, T1027I, V1176F)","count": len(data[data['lineage'] =='P.1'])},
    {"num":24,"name":"P.3 = P.3 (mutations:E484K, N501Y, D614G, P681H)","count":len(data[data['lineage'] =='P.3'])},
    {"num":25,"name":"S_GENE_DELETION = Variant virus with deletion in S-gene (defined by mutation: del 69-70 or by negative S-gene RT-PCR)","count":0},
    {"num":26,"name":"UNK = Sequence information unknown or not available","count":len(data[data['lineage'] =='None'])},
    {"num":27,"name":"VARIANT_OTHER = Novel variant of potential concern. Provide details in VirusVariantOther","count":0},
    {"num":28,"name":"WILD_TYPE = None of the variants described for this variable","count":0},
    {"num":29,"name":"Y453F = Y453F associated with farmed minks; defined by mutation: Y453F","count":0}
]

cur_year = datetime.date.today().isocalendar()[0]
cur_week = datetime.date.today().isocalendar()[1]

#Declare template variables
context = {
    'lineage_data_table': lineage_data_table,
    'cur_week': cur_week - 1,
    'cur_year': cur_year
    }

#Render automated report
template.render(context)
template.save(f'{output_path}/TESSY_{cur_week-1}ned_{cur_year}_NMRL.docx')
