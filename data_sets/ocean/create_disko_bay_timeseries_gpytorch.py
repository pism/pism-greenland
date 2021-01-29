#!/usr/bin/env python
# Copyright (C) 2020-21 Andy Aschwanden

from datetime import datetime
import gpytorch
import torch
import numpy as np
import pandas as pd
import pylab as plt
import time
import os

def toDecimalYear(date):
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year + 1, month=1, day=1)
    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed / yearDuration

    return date.year + fraction



if __name__ == "__main__":

    # depths to average over
    depth_min = 225
    depth_max = 275
    # depth for freezing point calculation
    depth = 250
    salinity = 34


    step = 1. / 12
    decimal_time = np.arange(1980, 2021 + step, step)

    mo_df = pd.read_csv("disko_bay_motyka.csv")

    omg_df = pd.read_csv("disko_bay_omg_axctd.csv", na_values=-99.0).dropna()
    omg_df = omg_df[(omg_df["Depth"] <= depth_max) & (omg_df["Depth"] >= depth_min)]
    omg_time = pd.to_datetime(omg_df.Date, format="%m/%d/%Y %H:%M:%S")
    omg_df["Date"] = omg_time
    omg_df = omg_df.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_df = omg_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    omg_time = [toDecimalYear(d) for d in omg_df.index]
    omg_df["Date"] = omg_time

    ices_df = pd.read_csv("disko_bay_ices.csv")
    ices_df = ices_df[(ices_df.Depth >= depth_min) & (ices_df.Depth <= depth_max)].reset_index(drop=True)
    ices_time = pd.to_datetime(ices_df.Date, format="%Y/%m/%d %H:%M:%S")
    ices_df["Date"] = ices_time
    ices_df = ices_df.set_index("Date")
    ices_df = ices_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    ices_time = [toDecimalYear(d) for d in ices_df.index]
    ices_df["Date"] = ices_time
    ices_df.to_csv("disko_bay_ices_depth_averaged.csv")

    holl_df = pd.read_csv("disko_bay_xctd_holland.csv", parse_dates=[0])
    holl_df = holl_df[(holl_df.Depth >= depth_min) & (holl_df.Depth <= depth_max)].reset_index(drop=True)
    holl_df = holl_df.set_index("Date")
    holl_df = holl_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    holl_time = [toDecimalYear(d) for d in holl_df.index]
    holl_df["Date"] = holl_time
    holl_df.to_csv("disko_bay_ices_depth_averaged.csv")


    X_mo = mo_df.Date.values.reshape(-1, 1)
    X_ices = ices_df.Date.values.reshape(-1, 1)
    X_omg = omg_df.Date.values.reshape(-1, 1)
    X_holl = holl_df.Date.values.reshape(-1, 1)

    y_mo = mo_df.Temperature.values
    y_ices = ices_df.Temperature.values
    y_omg = omg_df.Temperature.values
    y_holl = holl_df.Temperature.values

    all_df = pd.concat([mo_df, ices_df, omg_df, holl_df])
    all_df = all_df.sort_values(by="Date")
    
    X = all_df.Date.values.reshape(-1, 1)
    y = all_df.Temperature.values
    X_new = decimal_time[:, None]

    # We will use the simplest form of GP model, exact inference
    class ExactGPModel(gpytorch.models.ExactGP):
        def __init__(self, train_x, train_y, likelihood, cov):
            super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
            self.mean_module = gpytorch.means.ConstantMean()
            self.covar_module = gpytorch.kernels.ScaleKernel(cov)

        def forward(self, x):
            mean_x = self.mean_module(x)
            covar_x = self.covar_module(x)
            return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
        
    X_train = torch.tensor(X).to(torch.float)
    y_train = torch.tensor(np.squeeze(y)).to(torch.float)
    X_test = torch.tensor(X_new).to(torch.float)

    # initialize likelihood and model
    likelihood = gpytorch.likelihoods.GaussianLikelihood()

    cov = gpytorch.kernels.RBFKernel
    model = ExactGPModel(X_train, y_train, likelihood, cov())

    # Find optimal model hyperparameters
    model.train()
    likelihood.train()

    # Use the adam optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=0.1)  # Includes GaussianLikelihood parameters

    # "Loss" for GPs - the marginal log likelihood
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

    for i in range(500):
        # Zero gradients from previous iteration
        optimizer.zero_grad()
        # Output from model
        output = model(X_train)
        # Calc loss and backprop gradients
        loss = -mll(output, y_train)
        loss.backward()
        if i % 20 == 0:
            print(i, loss.item(), model.likelihood.noise.item())
        optimizer.step()


    # Get into evaluation (predictive posterior) mode
    model.eval()
    likelihood.eval()
    with torch.no_grad():  # , gpytorch.settings.fast_pred_var():
        # Draw n_samples
        n_samples = 10
        f_pred = model(X_test)
        samples = f_pred.sample(sample_shape=torch.Size([n_samples,]))
        
        # Initialize plot
        fig, ax = plt.subplots(1, 1)

        ax.plot(X_test.numpy(), samples.numpy().T, color='k', linewidth=0.5)

        # plot the data and the true latent function
        ax.plot(X_holl, y_holl, "o", color="#f4a582", ms=4, label="Observed (Holland)")
        ax.plot(X_mo, y_mo, "o", color="#b2182b", ms=4, label="Observed (Motyka)")
        ax.plot(X_ices, y_ices, "o", color="#92c5de", ms=4, label="Observed (ICES)")
        ax.plot(X_omg, y_omg, "o", color="#2166ac", ms=4, label="Observed (OMG)")

        ax.set_xlabel("Time")
        ax.set_ylabel("Temperature (Celsius)")
        ax.set_xlim(1980, 2021)
        ax.set_ylim(0, 5)
        plt.legend()
        fig.savefig(
            "disko-bay-temps.pdf"
        )

