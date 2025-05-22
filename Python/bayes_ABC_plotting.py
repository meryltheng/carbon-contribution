from bayes_ABC import *
import matplotlib
import random

def combine_data(num_samples,threshold_known, threshold_freq,max_it_number):

    folder_path = os.path.join(os.path.dirname(__file__))
    df_list = []

    for it_num in range(0,max_it_number+1):
        csv_file_name = os.path.join(folder_path, "bayes_samples", f"bayes_posteriors_thres{threshold_known}_{threshold_freq}_num{num_samples}_{it_num}.csv")

        if os.path.isfile(csv_file_name):
            output = pd.read_csv(csv_file_name)
            df_list.append(output)

    mega_DF = pd.concat(df_list)

    theta_RT_post_combined = mega_DF["theta_RT"].to_list()
    theta_MS_post_combined =  mega_DF["theta_MS"].to_list()
    theta_HS_post_combined =  mega_DF["theta_HS"].to_list()
    mean_post_combined = mega_DF["mean"].to_list()
    sd_post_combined = mega_DF["sd"].to_list()


    return  theta_RT_post_combined, theta_MS_post_combined, theta_HS_post_combined, mean_post_combined, sd_post_combined


def plot_threshold_densities(theta_RT_post, theta_MS_post, theta_HS_post, folder_path):
    density_RT = gaussian_kde(np.array(theta_RT_post))
    density_RT.covariance_factor = lambda : .25
    density_RT._compute_covariance()

    density_MS = gaussian_kde(np.array(theta_MS_post))
    density_MS.covariance_factor = lambda : .25
    density_MS._compute_covariance()

    density_HS = gaussian_kde(np.array(theta_HS_post))
    density_HS.covariance_factor = lambda : .25
    density_HS._compute_covariance()

    xs = np.linspace(0, 1, 200)
    fig, ax = plt.subplots(1,1, figsize=(10,5))

    cmap = matplotlib.colormaps['inferno']
    reversed_cmap = cmap.reversed()

    ax.fill_between(xs, density_RT(xs), alpha=0.2,color=reversed_cmap(1/5))
    ax.plot(xs, density_RT(xs),label="$\\theta_{RT}$",color=reversed_cmap(1/5))
    ax.fill_between(xs, density_MS(xs), alpha=0.2,color=reversed_cmap(2/5))
    ax.plot(xs,density_MS(xs),label="$\\theta_{MS}$",color=reversed_cmap(2/5))
    ax.fill_between(xs, density_HS(xs), alpha=0.2,color=reversed_cmap(3/5))
    ax.plot(xs, density_HS(xs),label="$\\theta_{HS}$",color=reversed_cmap(3/5))
    # ax.set_title("posterior thresholds")
    ax.set_xlabel("Threshold value")
    ax.set_ylabel("density")
    ax.legend()
    ax.set_xlim([0,1])

    ax.grid()

    plt.savefig( os.path.join(folder_path, "figures", "bayes_thresholds_posteriors.png"), bbox_inches='tight')
    plt.savefig( os.path.join(folder_path, "figures", "bayes_thresholds_posteriors.pdf"), bbox_inches='tight')

def plot_damage_normal_param_densities(mean_post,sd_post, folder_path):
    density_mean = gaussian_kde(np.array(mean_post))
    density_mean.covariance_factor = lambda : .25
    density_mean._compute_covariance()

    density_sd = gaussian_kde(np.array(sd_post))
    density_sd.covariance_factor = lambda : .25
    density_sd._compute_covariance()

    xs = np.linspace(-5, 6, 200)
    fig, ax = plt.subplots(1,1, figsize=(10,5))
    ax.fill_between(xs, density_mean(xs), alpha=0.1,color="#0066cc")
    ax.fill_between(xs, density_sd(xs), alpha=0.1,color="#00cc99")
    ax.plot(xs, density_mean(xs),label="$\mu$ mean",color="#0066cc")
    ax.plot(xs,density_sd(xs),label="$\sigma$ standard deviation",color="#00cc99")
    # ax.set_title("parameters for damage density normal")
    ax.set_xlabel("value")
    ax.set_ylabel("density")
    ax.legend()

    ax.grid()
    ax.set_xlim([-5,6])

    plt.savefig( os.path.join(folder_path, "figures", "bayes_damagenormal_posteriors.png"), bbox_inches='tight')
    plt.savefig( os.path.join(folder_path, "figures", "bayes_damagenormal_posteriors.pdf"), bbox_inches='tight')


def plot_parameter_set(index, theta_RT_post, theta_MS_post, theta_HS_post, mean_post,sd_post,folder_path ):
    xs = np.linspace(0, 1, 200)
    fig, ax = plt.subplots(1,1, figsize=(8,4))
    ax.plot(xs, truncated_norm_pdf(xs,mean_post[index], sd_post[index]),label="damage distribution",color="crimson",linewidth =4 )

    plt.axvline(x=0,ls='--',color="black")
    plt.axvline(x=1,ls='--',color="black")

    plt.axvline(x=theta_RT_post[index],ls='--',color="black")
    ax.annotate("$\\theta_{RT}$",
            xy=(theta_RT_post[index]-0.08, 0.15), xycoords='data',
            xytext=(1.5, 1.5), textcoords='offset points',
            fontsize=18)

    plt.axvline(x=theta_MS_post[index],ls='--',color="black")
    ax.annotate("$\\theta_{MS}$",
            xy=(theta_MS_post[index]-0.08, 0.15), xycoords='data',
            xytext=(1.5, 1.5), textcoords='offset points',
            fontsize=18)
    
    plt.axvline(x=theta_HS_post[index],ls='--',color="black") 
    ax.annotate("$\\theta_{HS}$",
            xy=(theta_HS_post[index]-0.08, 0.15), xycoords='data',
            xytext=(1.5, 1.5), textcoords='offset points',
            fontsize=18)

    cmap = matplotlib.colormaps['inferno']
    reversed_cmap = cmap.reversed()

    xs = np.linspace(0, theta_RT_post[index], 200)
    ax.fill_between(xs,truncated_norm_pdf(xs, mean_post[index], sd_post[index]),label="Relatively tolerant", alpha=0.3,color=reversed_cmap(1/5))

    xs = np.linspace(theta_RT_post[index],theta_MS_post[index], 200)
    ax.fill_between(xs,truncated_norm_pdf(xs, mean_post[index], sd_post[index]),label="Moderate susceptibility", alpha=0.3,color=reversed_cmap(2/5))

    xs = np.linspace(theta_MS_post[index],theta_HS_post[index], 200)
    ax.fill_between(xs,truncated_norm_pdf(xs, mean_post[index], sd_post[index]),label="High susceptibility", alpha=0.3,color=reversed_cmap(3/5))

    xs = np.linspace(theta_HS_post[index],1, 200)
    ax.fill_between(xs,truncated_norm_pdf(xs, mean_post[index], sd_post[index]),label="Extreme susceptibility", alpha=0.3,color=reversed_cmap(4/5))

    ax.legend(fontsize=16,framealpha=1)

    ax.grid()

    ax.tick_params(axis='both', which='major', labelsize=14)
    # ax.tick_params(axis='both', which='minor', labelsize=8)

    ax.set_xlabel("damage",fontsize=16)
    ax.set_ylabel("density",fontsize=16)

    plt.savefig( os.path.join(folder_path, "figures", f"bayes_sample_dis-{index}.png"), bbox_inches='tight')

def plot_posterior_across_included_samples(posterior,x_label,plot_title,xlims, folder_path, param_name):
    xs = np.linspace(xlims[0], xlims[1], 200)
    fig, ax = plt.subplots(1,1, figsize=(8,4))

    cmap = matplotlib.colormaps['inferno']
    reversed_cmap = cmap.reversed()

    step = int(np.floor(len(posterior)/10))
    for i in range(0,10):
        num = step*(i+1)
    
        density = gaussian_kde(np.array(posterior[0:num]))
        density.covariance_factor = lambda : .25
        density._compute_covariance()
        
        ax.plot(xs, density(xs), alpha=0.9, label=f"{num} samples",color=reversed_cmap((i+1)/10))

    ax.set_title(plot_title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("density")
    ax.legend()

    ax.grid()

    plt.savefig( os.path.join(folder_path , "figures", f"bayes_posterior_{param_name}.png"), bbox_inches='tight')


def plot_posterior_limited_samples(num_samples_included, posterior,x_label,plot_title,xlims, folder_path, param_name):
    xs = np.linspace(xlims[0], xlims[1], 200)
    fig, ax = plt.subplots(1,1, figsize=(8,4))

    cmap = matplotlib.colormaps['winter']
    reversed_cmap = cmap.reversed()
    for i in range(0,10):
        limited_posterior = random.sample(posterior,num_samples_included)
        density = gaussian_kde(np.array(limited_posterior))
        density.covariance_factor = lambda : .25
        density._compute_covariance()
        
        ax.plot(xs, density(xs), alpha=0.9, color=reversed_cmap((i+1)/10))

    ax.set_title(plot_title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("density")
    # ax.legend()

    ax.grid()

    plt.savefig( os.path.join(folder_path , "figures", f"bayes_posterior_{param_name}_{num_samples_included}.png"), bbox_inches='tight')


def plotting(num_samples,threshold_known,threshold_freq,max_it_number = 10, folder_path = os.path.join(os.path.dirname(__file__))):

    theta_RT_post, theta_MS_post, theta_HS_post, mean_post,sd_post = combine_data(num_samples,threshold_known, threshold_freq,max_it_number)

    print(len(theta_RT_post))

    plot_threshold_densities(theta_RT_post, theta_MS_post, theta_HS_post, folder_path)

    plot_damage_normal_param_densities(mean_post,sd_post, folder_path)

    return

    for index in [2,10]: #  [0,10,20,30,40,50]:  # range(0,10)  # 
        plot_parameter_set(index, theta_RT_post, theta_MS_post, theta_HS_post, mean_post,sd_post,folder_path )

    

    plot_posterior_across_included_samples(theta_RT_post,"RT threshold","Posterior of $\\theta_{RT}$ threshold",[0,1], folder_path, "thresholdRT")

    plot_posterior_across_included_samples(theta_MS_post,"MS threshold","Posterior of $\\theta_{MS}$ threshold",[0,1], folder_path, "thresholdMS")

    plot_posterior_across_included_samples(theta_HS_post,"HS threshold","Posterior of $\\theta_{HS}$ threshold",[0,1], folder_path, "thresholdHS")

    plot_posterior_across_included_samples(mean_post,"damage norm mean","Posterior of $\\mu$",[-5,2], folder_path, "mean")

    plot_posterior_across_included_samples(sd_post,"damage norm sd","Posterior of $\\sigma$",[0,10], folder_path, "standarddeviation")

    plot_posterior_limited_samples(25, mean_post,"damage norm mean","Posterior of $\\mu$ - 25 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(50, mean_post,"damage norm mean","Posterior of $\\mu$ - 50 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(100, mean_post,"damage norm mean","Posterior of $\\mu$ - 100 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(200, mean_post,"damage norm mean","Posterior of $\\mu$ - 200 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(500, mean_post,"damage norm mean","Posterior of $\\mu$ - 500 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(750, mean_post,"damage norm mean","Posterior of $\\mu$ - 750 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(1000, mean_post,"damage norm mean","Posterior of $\\mu$ - 1000 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(2000, mean_post,"damage norm mean","Posterior of $\\mu$ - 2000 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(5000, mean_post,"damage norm mean","Posterior of $\\mu$ - 5000 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(7500, mean_post,"damage norm mean","Posterior of $\\mu$ - 7500 random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(10000, mean_post,"damage norm mean","Posterior of $\\mu$ - 10k random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(12000, mean_post,"damage norm mean","Posterior of $\\mu$ - 12k random samples",[-5,2], folder_path, "mean")
    plot_posterior_limited_samples(15000, mean_post,"damage norm mean","Posterior of $\\mu$ - 15k random samples",[-5,2], folder_path, "mean")


if __name__ == "__main__":
    plotting(num_samples=10,threshold_known=0.5, threshold_freq=0.1,max_it_number = 2000)