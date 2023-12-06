#!/usr/bin/env python3
#-*- coding: utf-8 -*-

'''
-------------------------------------------------------------------------------
Copyright© 2021 Jules Devaux, Alice Harford, Tony Hickey. All Rights Reserved
Open Source script under Apache License 2.0
-------------------------------------------------------------------------------
'''

import os
current_path = os.path.dirname(os.path.abspath(__file__))


#========= PARAMETERS =========
TEMPERATURES=[18, 20, 24, 25, 26, 27, 30]# 
CHAMBER_VOLUME = 2 #ml
WINDOW = 10 #Moving averages
MAX_kPA = 24 #kPa
GRAPHING = False
_saving = True # Save summary file and stats
FOLDER = f'{current_path}/CSV/Leak_OXPHOS_mtMP_JO2_PO2_CSV_Original' # Folder containing all csvs
SUBFOLDER_NAME='Leak_OXPHOS_mtMP_JO2_PO2_'

#=============================


import sys, pickle, datetime
import numpy as np
import pandas as pd
import pingouin as pg
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


def hill(S, Vmax, Km, n=1):
	S=np.abs(S)
	return (Vmax*S**n)/(Km+S**n)


def fit_hill_curve(x, y):
	'''
	x_y is df with x as index, y as values
	'''
	Vmax = y.max()
	xmax, xmin = x[-1], x[0]
	dx = x[-1]-x[0]

	# parameters as [max, min, dt, ]
	param_initial = [Vmax, 0.5*dx, 1]
	param_bounds = ([Vmax-0.5*Vmax, xmin*0.1, 0.01],
					[Vmax+0.5*Vmax, xmax*10, 100])
	bounds = [(Vmax-0.5*Vmax, Vmax+0.5*Vmax),
			(xmin, xmax),
			(0.01, 100)]


	popt, _ = curve_fit(
		hill, x, y,
		p0 = param_initial,
		bounds=param_bounds)

	Vmax, pP50, _hill =popt[0], popt[1], popt[2]
	y_fit = hill(x, Vmax, pP50, _hill)


	# Goodness of fit:
	# residual sum of squares
	ss_res = np.sum((y - y_fit) ** 2)
	ss_tot = np.sum((y - np.mean(y)) ** 2)
	r2 = 1 - (ss_res / ss_tot)

	predicted=pd.Series(y_fit, name='predicted')

	return Vmax, pP50, _hill, predicted, r2


def get_JO2(PO2, mass=2, Kpa_to_pmol=True):
	'''
	Recalculates JO2 from PO2.
	Expects PO2 as series and in kPa.
	Requires weight in mg (float)
	'''

	# Moving average
	PO2=PO2.rolling(window=int(WINDOW/4)).mean()

	# Calculate JO2
	# Expressed in kPa/s/mg
	JO2=((((PO2.shift())-PO2)*CHAMBER_VOLUME)/(mass*2))*1000

	# Exponential moving average to smooth it all.
	JO2=JO2.ewm(span=WINDOW, adjust=False).mean()

	# Convert to pmol/s/mg
	if Kpa_to_pmol is True:
		JO2=JO2*(224.47/20.37)

	return JO2




def extract_csvs():
	# All chambers saved as dict
	chambers = []

	for temperature in TEMPERATURES:

		# List all csv in temperature folders
		_path = f"{FOLDER}/{SUBFOLDER_NAME}{temperature}/"
		csvs=[file for file in os.listdir(_path) if '.csv' in file]

		for csv in csvs:
			print(f'Exctracting {csv}...')
			df=pd.read_csv(f"{_path}{csv}", low_memory=False, encoding= 'unicode_escape')
			
			# Split each chamber:
			for chb in ['A', 'B']:
				try:
					# Retrieve weight
					WGH_col=[c for c in df.columns if f"Weight {chb}" in c][0]
					weight=float(df.loc[:,WGH_col].iloc[0])

					# Retrieve and select columns for that specific chamber
					PO2_col=[c for c in df.columns if f"{chb}: O2 pres" in c][0]
					JO2_col=[c for c in df.columns if f"{chb}: O2 flux" in c][0]
					DYm_col=[c for c in df.columns if f"{chb}: Amp [" in c][0]

					# Create DataFrame for that chamber
					cdf=df.loc[:,[PO2_col,JO2_col,DYm_col]].dropna()
					cdf=cdf.rename(columns={PO2_col:'PO2',
											JO2_col:'JO2',
											DYm_col:'DYm'})


					# Recalculate JO2 from PO2
					cdf['JO2_calc']=get_JO2(cdf['PO2'], weight)

					cdf['PO2']=cdf['PO2']-cdf['PO2'].min()
					cdf['JO2']=cdf['JO2']-cdf['JO2'].min()
					cdf['JO2_calc']=cdf['JO2_calc']-cdf['JO2_calc'].min()


					# Define OXPHOS - LEAK
					if chb[0] == 'A':
						_state = 'OXPHOS'
					else: _state = 'LEAK'

					chambers.append({
						'filename':csv,
						'temperature':temperature,
						'chamber':chb,
						'state':_state,
						'df':cdf,
						'mass':weight
						})


				except Exception as e:print(f"ERROR: {csv}: chamber {chb}: {e}")
			
	return chambers


def fundamentals(chambers, _saving=True):

	_dicts=[] # Main list to create the final df
	for chamber in chambers:

		_dicts.append({
			'filename':f"{chamber['filename']}_{chamber['state']}",
			'PO2_max':chamber['df']['PO2'].max(),
			'PO2_min':chamber['df']['PO2'].min(),
			'JO2_max':chamber['df']['JO2'].max(),
			'JO2_min':chamber['df']['JO2'].min(),
			'dt':len(chamber['df'])
		})

	df=pd.DataFrame(_dicts)
	summary=df.agg(['mean','sem'])
	if _saving:
		df.to_csv('fundamentals.csv')
		print('saved fundamentals')

	return df


def process_all_chambers(chambers):
	'''
	index as time, columns = JO2, PO2, DYm
	'''

		
	_dicts=[] #temprary list to avoid overwriting chambers list
	for chamber in chambers:
		try:

			print(f"Processing {chamber['filename']}")

			

			# GRAPH
			#plt.scatter(X, predicted, label='pred')
			# plt.scatter(chamber['df'].index, chamber['df']['DYm'].values, label='expe')
			# plt.title(chamber['filename'])
			# plt.show()

			# ============= JO2 PARAMETERS AND FITTED HILL CURVE ============
			# Observed P50 from raw JO2
			try:
				JO2_max=chamber['df']['JO2'].max()
				JO2_50=JO2_max/2
				_portion=chamber['df'].loc[chamber['df']['JO2']<(JO2_50)]
				oP50=_portion['PO2'].iloc[0]
			except Exception as e:
				print(f"ERROR: observed_P50: {e}")
				oP50=float('nan')

			# need to have df with PO2 index and JO2 column
			tdf=chamber['df'].loc[:,['PO2','JO2']]
			tdf=tdf.drop_duplicates(subset='PO2', keep='first').set_index('PO2').sort_index()

			# Fit curve to experimental PO2 
			fJO2max, fP50, _hill, fitted, r2 = fit_hill_curve(tdf.index.values, tdf['JO2'].values)

			# Extend curve to 0 -> 100% PO2
			# Standardise PO2 for comaprison between samples
			X=np.arange(0, MAX_kPA, 0.01)
			predicted=hill(X, fJO2max, fP50, _hill)
			predicted=pd.DataFrame(predicted, index=X, columns=['JO2'])

			# Calculate area under the curve from Standardised
			auc=np.trapz(predicted['JO2'], predicted.index, 0.1)

			# JO2max at 100% saturation 
			kPa21_JO2max=predicted['JO2'].max()

			# Extract P50 from fitted JO2
			try:
				JO2_50=kPa21_JO2max/2
				_portion=predicted.loc[predicted['JO2']<JO2_50]
				cP50=_portion.index.values[-1]
			except Exception as e:
				print(f"ERROR: calculated_P50: {e}")
				cP50=float('nan')



			# ============= MEMBRANE POTENTIAL ==========
			dym=chamber['df'].loc[:,['PO2','DYm']]
			dym=dym.drop_duplicates(subset='PO2', keep='first'
					).set_index('PO2').sort_index()
			dym=dym.loc[(dym.index<21)]
			dym=(2-dym)/chamber['mass']
			dym=dym.rolling(window=WINDOW).mean().dropna()


			# Fit curve to experimental PO2 
			fDYmmax, fP50_DYm, DYm_hill, DYm_fitted, DYm_r2 = fit_hill_curve(dym.index.values, dym['DYm'].values)

			# Extend curve to 0 -> 100% PO2
			# Standardise PO2 for comaprison between samples
			DYm_predicted=hill(X, fDYmmax, fP50_DYm, DYm_hill)
			DYm_predicted=pd.DataFrame(DYm_predicted, index=X, columns=['DYm'])

			# Calculate area under the curve from Standardised
			DYm_auc=np.trapz(DYm_predicted['DYm'], DYm_predicted.index, 0.1)

			# DYmmax at 100% saturation 
			kPa21_DYmmax=DYm_predicted['DYm'].max()

			# Extract P50 from fitted DYm
			try:
				DYm_50=kPa21_DYmmax/2
				_portion=DYm_predicted.loc[DYm_predicted['DYm']<DYm_50]
				cP50_DYm=_portion.index.values[-1]
			except Exception as e:
				print(f"ERROR: DYm_calculated_P50: {e}")
				cP50_DYm=float('nan')



			# Update parameters for chamber
			chamber.update(
				{'JO2_observed_P50':oP50,
				'JO2max_fitted': fJO2max,
				'JO2_fitted_P50':fP50,
				'JO2_fitted_hill':_hill,
				'JO2max_at_21kPa':kPa21_JO2max,
				'JO2_calc_P50':cP50,
				'JO2 goodness of fit':r2,
				'JO2 area under curve':auc,
				'JO2 predicted':predicted,
				'DYmmax_fitted': fDYmmax,
				'DYm_fitted_P50':fP50_DYm,
				'DYm_fitted_hill':DYm_hill,
				'DYmmax_at_21kPa':kPa21_DYmmax,
				'DYm_calc_P50':cP50_DYm,
				'DYm goodness of fit':DYm_r2,
				'DYm area under curve':DYm_auc,
				'DYm predicted':DYm_predicted,
				})

			_dicts.append(chamber)

			# GRAPH
			plt.scatter(DYm_predicted.index, DYm_predicted['DYm'], label='pred', c='blue')
			plt.scatter(dym.index, dym.values, label='expe', c='red')
			plt.title(f"{chamber['filename']} - {chamber['state']}")
			plt.show()

		except Exception as e:print(e)

	return _dicts


def create_summary(chambers):

	# Check outliers
	outliers=[]
	for chamber in chambers:
		if (float(chamber['JO2 goodness of fit'])<=0.5) or (float(chamber['DYm goodness of fit'])<=0.5):
			outliers.append(chamber['filename'])
			# outliers.append({
			# 	'filename':chamber['filename'],
			# 	'chamber': chamber['chamber']
			# 	})

		print("============= OUTLIERS ==============")
		for o in outliers:
			print(o)
			#print(f"{o['filename']} - chamber {o['chamber']}")
		print("\n\n")

		# Create summary table
		summary=[]
		for i in range(len(chambers)):
			if chambers[i]['filename'] not in outliers:
				row={k:v for k,v in chambers[i].items() if k not in ['df','PO2','JO2 predicted','JO2','JO2_calc','DYm', 'DYm predicted']}
				row=pd.DataFrame(row, index=[i])
				summary.append(row)
		summary=pd.concat(summary)

		# Save summary
		if _saving is True:
			summary.to_csv('summary.csv')

		return summary


def do_stats(df):
	
	# Create Excel template
	_now=datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
	fileName=f'Summary_{_now}.xlsx'
	writer=pd.ExcelWriter(fileName, engine='xlsxwriter')

	# First, add data to first tab
	df.to_excel(writer, sheet_name='main') 


	between=['temperature','state']
	parameters=[c for c in df.columns if df[c].dtypes != 'object' and c not in between]


	#===== Do Averages and SEM
	summary=df.groupby(between).agg(['mean','sem'])
	summary.to_excel(writer, sheet_name='averages') 


	#===== Do two-Way ANOVA ====
	aov_row=0 #Keep track of rows to append ANOVA tests
	hoc_row=0

	# Loop through parameter for ANOVA
	for parameter in parameters:
		aov=pg.anova(df, dv=parameter,
				between=between,
				detailed=True)
		aov['Parameter']=parameter

		# Append to Excel, spacing of 3 rows
		aov.to_excel(writer, sheet_name="ANOVAs", startrow=aov_row , startcol=0)   
		aov_row+=(len(aov.index)+3)

		# ==== POST-HOCs
		groups=[{'groups': 'state', 'between': 'temperature'},
				{'groups': 'temperature', 'between': 'state'}]
		for g in groups:
			for group, gdf in df.groupby(g['groups']):
				# Do the test
				ph=pg.pairwise_tukey(data=gdf, dv=parameter, between=g['between'])
				ph['Parameter']=parameter
				ph['group']=group
		
				# Append to Excel
				ph.to_excel(writer, sheet_name="PostHocs", startrow=hoc_row , startcol=0)   
				hoc_row+=(len(ph.index)+3)


	writer.save()
	print(f'Saved {fileName}')



def main():
	chambers=extract_csvs()
	_fundamentals=fundamentals(chambers)
	chambers=process_all_chambers(chambers)
	with open('chambers.pkl','wb') as f:
		pickle.dump(chambers,f)
	df=create_summary(chambers)
	df=pd.read_csv('summary.csv', header=0).drop(columns='Unnamed: 0')
	do_stats(df)


if __name__=='__main__':
	main()
