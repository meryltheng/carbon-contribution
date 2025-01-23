# TODO: to generate csv with damages selected, so that Meryl can use them to run things.
# Need to think whether to sample those with known damages or not.
# probably don't sample -- consider if I had to use another data source, which may also have ratings...


from bayes_ABC_plotting import *


def estimate_damage(score, theta_RT, theta_MS, theta_HS,mean, sd):
    lowerbound = 1
    upperbound = 0
    if "ES" in score:
        upperbound = 1
    elif "HS" in score:
        upperbound = theta_HS 
    elif "MS" in score:
        upperbound = theta_MS 
    elif "RT" in score:
        upperbound = theta_RT
    else:
        print("something wrong about the score")
        exit(1)
    
    if "RT" in score:
        lowerbound = 0
    elif "MS" in score:
        lowerbound = theta_RT
    elif "HS" in score:
        lowerbound = theta_MS
    elif "ES" in score:
        lowerbound = theta_HS 
    else:
        print("something wrong about the score")
        exit(1)


    # sample damage from truncated normal distribution
    damage_est = truncated_norm(mean, sd,lowerbound,upperbound)
    
    return  damage_est

def damage_sample(row):
    if not pd.isna(row['Rating']) and pd.isna(row['Damage']):
        row['Damage'] = estimate_damage(row['Rating'], theta_RT_post[index], theta_MS_post[index], theta_HS_post[index],mean_post[index], sd_post[index])
    else:
        pass
    return row['Damage']


num_samples=10
threshold_known=0.5
threshold_freq=0.1
max_it_number = 2000


theta_RT_post, theta_MS_post, theta_HS_post, mean_post,sd_post = combine_data(num_samples,threshold_known,threshold_freq,max_it_number)

folder_path = os.path.join(os.path.dirname(__file__))

num_damage_samples = 1000

for i in range(1,num_damage_samples+1):

    df_combined = pd.read_csv( os.path.dirname(__file__) + '/MR_impact_base.csv',encoding='mbcs')

    index = np.random.randint(0, high=len(theta_RT_post))
    print(index)

    df_combined['Damage'] = df_combined.apply(damage_sample, axis=1)

    df_combined = df_combined.drop('Notes', axis=1)
    df_combined = df_combined.drop('Reference', axis=1)
    df_combined = df_combined.drop('Rating', axis=1)

    df_combined.to_csv(os.path.join(folder_path, "MR_samples", f"MR_impact_{i}.csv"), index=False)
