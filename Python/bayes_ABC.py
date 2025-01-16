import numpy as np
from scipy.stats import truncnorm
import pandas as pd
import os 
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import rc, cm
rc('text', usetex=True)
rc('font', **{'family': 'sans-serif'})
from scipy.stats import gaussian_kde 
from scipy.stats import norm 
import pickle
import csv
import argparse

def simulated_score_given_damage(damage,theta_RT, theta_MS, theta_HS):
    if damage< theta_RT:
        return "RT"
    elif damage< theta_MS:
        return "MS"
    elif damage < theta_HS:
        return "HS"
    else:
        return "ES"

def truncated_norm(mean, sd, lower=0,upper=1):
    a_trunc  = lower  # truncation
    b_trunc = upper
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html
    loc = mean
    scale = sd
    a, b = (a_trunc - loc) / scale, (b_trunc - loc) / scale # the way a,b are defined by scipy
    return truncnorm.rvs(a, b,loc=loc,scale=scale, size=1)[0] # sample one

def truncated_norm_pdf(x, mean, sd, lower=0,upper=1):
    a_trunc  = lower  # truncation
    b_trunc = upper
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html
    loc = mean
    scale = sd
    a, b = (a_trunc - loc) / scale, (b_trunc - loc) / scale # the way a,b are defined by scipy
    return truncnorm.pdf(x, a, b,loc=loc,scale=scale)

def simulated_score(theta_RT, theta_MS, theta_HS,mean, sd):
    # sample damage from truncated normal distribution
    damage_est = truncated_norm(mean, sd)
    
    return  simulated_score_given_damage(damage_est,theta_RT, theta_MS, theta_HS)

def score_difference(first_score,second_score):
    numerical_conversion = {'RT':0, 'MS':1, 'HS':2, 'ES':3, 'RT-MS-HS':1, 'MS-HS':1.5, 'HS-ES':2.5,'RT-MS-HS-ES':1.5, 'MS-HS':1.5, 'MS-HS-ES':2, 'RT-MS':0.5}

    first = numerical_conversion[first_score]
    second = numerical_conversion[second_score]
    return abs(first-second)


def damage_score_difference(actual_damage,actual_score,theta_RT, theta_MS, theta_HS):
    if actual_score=="RT" and actual_damage > theta_RT:
        return actual_damage-theta_RT 
    if actual_score=="MS" and actual_damage>theta_MS:
        return actual_damage-theta_MS 
    if actual_score=="HS" and actual_damage>theta_HS:
        return actual_damage-theta_HS
    if actual_score =="ES" and actual_damage < theta_HS:
        return theta_HS - actual_damage
    if actual_score =="HS" and actual_damage < theta_MS:
        return theta_MS - actual_damage
    if actual_score =="MS" and actual_damage < theta_RT:
        return theta_RT - actual_damage
    if actual_score=="RT-MS": # this means that the predicted score should have been above MS
        return actual_damage-theta_MS
    if actual_score =="HS-ES": # this means that the predicted score was below HS
        return theta_MS-actual_damage
    else:
        print(actual_score)
        exit(1)

def difference_full_data(known_scores, known_damages,theta_RT, theta_MS, theta_HS):
    difference = 0
    for score,damage in zip(known_scores,known_damages):
        pred_score= simulated_score_given_damage(damage,theta_RT, theta_MS, theta_HS)
        if pred_score in score:
            pass # difference of zero
        else:
            # difference+= score_difference(score, pred_score)

            # rather than the score conversion, the distance should be more like some kind of distance to the threshold that would make it correct? 
            difference += damage_score_difference(damage,score,theta_RT, theta_MS, theta_HS)

    return difference


def difference_distribution(known_freq,num_samples,theta_RT, theta_MS, theta_HS,mean, sd):
    list_of_simulated_scores = [simulated_score(theta_RT, theta_MS, theta_HS,mean, sd) for i in range(num_samples)]

    ratings_list = ['RT','MS','HS','ES']

    simulated_freq = {'RT':0, 'MS':0, 'HS':0,'ES':0}
    for val in list_of_simulated_scores:
        for sc in ratings_list:
            if sc in val:
                simulated_freq[sc] +=1
    

    simulated_freq_density = {sc:simulated_freq[sc]/num_samples  for sc in ratings_list}

    # https://en.wikipedia.org/wiki/Statistical_distance 
    # something like statistical distance
    difference = 0
    for sc in ['RT','MS','HS','ES']:
        difference+=abs(known_freq[sc]-simulated_freq_density[sc])

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wasserstein_distance.html


    return difference/2



def ABC(num_samples, threshold):

    # posterior distributions
    theta_RT_post = []
    theta_MS_post = []
    theta_HS_post = []
    mean_post = []
    sd_post = []

    # getting the data
    folder_path = os.path.join(os.path.dirname(__file__))
    df = pd.read_csv( os.path.dirname(__file__) + '/data_myrtle_rust_susceptibility.csv')
    not_null_mask = df.notnull().all(axis=1)
    df_OG = df[not_null_mask]
    df_new = df[~not_null_mask]

    # data set 1
    x_name = "Rating"
    y_name = "Damage"
    x = df_OG[x_name].tolist()
    y = df_OG[y_name].tolist()

    # data set two:
    x2 = df[x_name].tolist()
    ratings_list = ['RT','MS','HS','ES']
    freq = {'RT':0, 'MS':0, 'HS':0,'ES':0}
    total_samples = 0
    for val in x2:
        for sc in ratings_list:
            if sc in val:
                freq[sc] +=1

    total_samples = sum([freq[sc] for sc in ratings_list])

    freq_density = {sc:freq[sc]/total_samples  for sc in ratings_list}
    num_tries = 0

    for i in range (0, num_samples):
        distance = threshold+1
        while distance > threshold: # while distance is too big, sample:
            num_tries+=1
            # sample from priors
            thresholds_ordered = False 
            while not thresholds_ordered:
                theta_RT = truncated_norm(mean=1/4, sd=1)
                theta_MS = truncated_norm(mean=2/4, sd=1)
                theta_HS = truncated_norm(mean=3/4, sd=1)
                if theta_RT<0 or theta_MS<0 or theta_HS<0:
                    exit(1)
                if theta_RT < theta_MS and theta_MS<theta_HS:
                    thresholds_ordered = True

            # mean and sd should probably from an inverse gamma 
            mean = np.random.normal(loc=0,scale=2)
            sd = np.random.uniform(0,10)

            # generate simulated data and compare with known data
            # compare with combination data (score + known damage)
            difference_known = difference_full_data(x, y,theta_RT, theta_MS, theta_HS)
            # TODO there is only 11 of these, so...actually, it would be good if like at least half of them were right??? 

            # plus compare with score only damage - frequency distribution
            difference_freq = difference_distribution(freq_density,total_samples,theta_RT, theta_MS, theta_HS,mean, sd)

            # sum up
            distance = difference_known+difference_freq

            # print(f"{difference_known}, {difference_freq}, {distance}")

        theta_RT_post.append(theta_RT)
        theta_MS_post.append(theta_MS)
        theta_HS_post.append(theta_HS)
        mean_post.append(mean)
        sd_post.append(sd)

        

        print(f"{difference_known}, {difference_freq}, {distance}")

    acceptance_rate = num_samples/num_tries
    print(f"Acceptance rate: {acceptance_rate}")



    return theta_RT_post, theta_MS_post, theta_HS_post, mean_post,sd_post

def one_iteration(folder_path, num_samples, threshold, it_num):
    csv_file_name = os.path.join(folder_path, "bayes_samples", f"bayes_posteriors_thres{threshold}_num{num_samples}_{it_num}.csv")

    if os.path.isfile(csv_file_name):
        return # skip this iteration if file already computed and exists

    theta_RT_post, theta_MS_post, theta_HS_post, mean_post,sd_post = ABC(num_samples, threshold)

    with open(csv_file_name, 'w',newline='') as f:

        wr = csv.writer(f)
        wr.writerow(["theta_RT","theta_MS", "theta_HS", "mean", "sd"])
        for i in range(num_samples):
            wr.writerow([theta_RT_post[i], theta_MS_post[i], theta_HS_post[i], mean_post[i],sd_post[i]])


def main(folder_path = os.path.join(os.path.dirname(__file__))):

    num_samples = 10 # better to run in smaller batches, which will get done in a reasonable time, and then combine them all    

    threshold =0.7

    for it_num in range(0,10):
        one_iteration(folder_path, num_samples, threshold, it_num)
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--it_num', '-num', default="")
    args = parser.parse_args()
    it_num = args.it_num
    if it_num == "":  # i.e., default:
        main()
    else:
        folder_path = os.path.join(os.path.dirname(__file__))
        num_samples = 10
        threshold =0.7
        one_iteration(folder_path, num_samples, threshold, int(it_num))