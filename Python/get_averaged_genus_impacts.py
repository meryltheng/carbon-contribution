import os 
import pandas as pd
import numpy as np

num_samples=10
threshold_known=0.5
threshold_freq=0.1
max_it_number = 2000

folder_path = os.path.join(os.path.dirname(__file__))

num_damage_samples = 1000

combined_list = {"Eucalyptus":[],"Syzygium":[],"Rhodamnia":[]}


for i in range(1,num_damage_samples+1):
    df = pd.read_csv( os.path.join(os.path.dirname(__file__), "MR_samples",  f"MR_impact_{i}.csv"),encoding='mbcs')
    df = df.drop('Species', axis=1)

    df = df.groupby(['Genus'])[['Damage']].agg('mean').reset_index()

    for genus in combined_list.keys():
        result = df.loc[df['Genus'] == genus]
        result = result['Damage'].tolist()
        combined_list[genus].append(result[0])

medians = {}
quantile025 = {}
quantile975 = {}
for genus in combined_list.keys():
    medians[genus] = np.median(combined_list[genus])
    quantile025[genus] = np.quantile(combined_list[genus],0.025)
    quantile975[genus] = np.quantile(combined_list[genus],0.975)

print(medians)
print(quantile025)
print(quantile975)

genus_list = list(set(df.loc[:,"Genus"]))
genus_list = [x.capitalize() for x in genus_list ]
genus_list = sorted(list(set(genus_list)))
complete_list = {x:[] for x in genus_list}

print(genus_list)

for i in range(1,num_damage_samples+1):
    df = pd.read_csv( os.path.join(os.path.dirname(__file__), "MR_samples",  f"MR_impact_{i}.csv"),encoding='mbcs')
    df = df.drop('Species', axis=1)

    df['Genus'] = df['Genus'].str.capitalize()

    df = df.groupby(['Genus'])[['Damage']].agg('mean').reset_index()

    for genus in genus_list:
        result = df.loc[df['Genus'] == genus]
        result = result['Damage'].tolist()
        complete_list[genus].append(result[0])


medians = {}
quantile025 = {}
quantile975 = {}
print("Genus & Median &  0.025 quantile &  0.975 quantile \\\\")
for genus in genus_list:
    medians[genus] = np.median(complete_list[genus])*100
    quantile025[genus] = np.quantile(complete_list[genus],0.025)*100
    quantile975[genus] = np.quantile(complete_list[genus],0.975)*100
    print(f"{genus} & {medians[genus]:.1f}\%  & {quantile025[genus]:.1f}\% &  {quantile975[genus]:.1f}\% \\\\ ")

# data = {'Genus': genus_list,
#         'Median': [medians[x] for x in genus_list],
#         '0.025 quantile':[quantile025[x] for x in genus_list],
#          '0.975 quantile': [quantile975[x] for x in genus_list] }

# df = pd.DataFrame(data)

# print(df.to_latex(index=False,

#                   float_format="{:.1f}\%".format,

# ))