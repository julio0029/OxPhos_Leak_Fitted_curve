#!/usr/bin/env python3
#-*- coding: utf-8 -*-
'''
-------------------------------------------------------------------------------
Copyright© 2021 Jules Devaux, Alice Harford, Tony Hickey. All Rights Reserved
Open Source script under Apache License 2.0
-------------------------------------------------------------------------------
'''


import os
current_path = f"{os.path.dirname(os.path.abspath(__file__))}/"



# ======  PARAMETERS  =====================

TEMPERATURES=[24,26,28,30,32,34]
DATAPATH=f"{current_path}"
REPO_PATH="/Users/julesdevaux/Dropbox/Alice Harford/PhD/Wrasse_Numbers_Dimensions/Wrasse_Dimension_Measurements_2020-2022.xlsx"
REPO_TABNAME = "Parrot_Wrasse_Chamber_Tisue"

ATP_calib = {'PMG':0, #mM
			'ATP_1':2.5, 
			'ATP_2':5}

ADP_calib = {'PMG':0, #mM
			'ADP_1':2.5,
			'ADP_2':5}

MgCl2_calib = {'MgG': 0, #mM
			'MgCl2_1':62.5,
			'MgCl2_2':125}

MAX_kPA=21


OUTLIERS=['ATP_BandedWrasse_Trial2_Ca2+_20oCWonderwoman_25.5.22.csv',	# Amp saturated in chamber B
		'ATP_BandedWrasse_Trial7_Ca2+_20oCMystique_8.6.22.csv',			# Amp saturated in both chambers
		'ATP_BandedWrasse_Trial6_Ca2+_25oCMystique_6.6.22.csv',			# Amp saturated in chamber A
		'ATP_BandedWrasse_Trial7_Ca2+_25oCWonderwoman_8.6.22.csv',		# Amp saturated in both chambers
		'ATP_BandedWrasse_Trial6_Ca2+_26oCMystique_12.6.22.csv',		# Amp saturated in chamber B
		'ATP_BandedWrasse_Trial5_Ca2+_27oCMystique_12.6.22.csv',		# Amp saturated in both chambers
		'ATP_BandedWrasse_Trial2_Ca2+_27oCMystique_25.5.22.csv',		# Amp saturated in both chambers
		'ATP_BandedWrasse_Trial7_Ca2+_30oCSuperman_8.6.22.csv',			# Amp saturated in both chambers
		'ATP_BandedWrasse_Trial6_Ca2+_30oCSuperman_6.6.22.csv',			# Amp saturated in both chambers for MgCl2 calibration
		'PW27_28oC_Trial4_Superman_HeronIsland_29.11.22.csv',			# MgCl2 calib not linear
		'PW29_28oC_Trial5_Superman_HeronIsland_29.11.22.csv',			# MgCl2 calib not linear
		]

#==========================================


'''
Code steps:

1) Create chamber list
	- Load repository for filename, conditions, mass
	- Loop through all csvs
		- creat chamber dict -> chamber class

2) Exctract ATP and ADP standard curves
	- sc for each file
	- as well as create df for each temperature and extract sc

3) Calibrate Amp to ATP

4) Extract all stable states

5) Create ATP profiles with anoxia
	- JO2 and ATP rates on PO2

6) Stats


'''


# Constants:
F = 96485.33212 #C mol−1
R = 8.314472 #J.K-1
Z = 1 #Valence of TMRM / safranine 
#=========================================================================================



#////////////////////////////////////// SCRIPT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
import time
import pandas as pd
import numpy as np
import pickle
from matplotlib import pyplot as plt
import pingouin as pg
from sklearn import linear_model
from scipy.optimize import curve_fit
import statsmodels.api as sm


FLUORESCENCE={
	'ATP': ['ATP', 
		'ADP',
		'MgG', 'MAGNESIUM GREEN'],

	'DYm': ['DYm',
		'Safr', 'SAFRANIN',
		'TMRM', 'TMRE'],

	'ROS': ['ROS',
		'AuR', 'AmplexRed'],

	'NADH': ['NADH', 'NAD'],

	'CytC': ['CytC', 'Cytochrome c'],

	'FADH': ['FADH', 'FAD'],

		}



def linear_fit(x, y, steps=0.1):

	#Get standard curve
	lm=linear_model.LinearRegression()
	_x=np.array(x).reshape(-1,1)

	model=lm.fit(_x,y)
	r2=lm.score(_x,y)
	slope=lm.coef_[0]
	intercept=lm.intercept_

	X=np.arange(min(x), max(x), steps)
	Y=model.predict(X.reshape(-1,1))
	
	predicted=pd.DataFrame(Y, index=X)

	return float(slope), float(intercept), r2, predicted


class ATP():
	def __init__(self, df=None, **kwargs):
		if df is not None:
			self.df=df

		self.__defaultParams__()

		for k, v in kwargs.items():
			setattr(self, k, v)


	def __defaultParams__(self):
		self.adtp='ADP'
		self.ATP_calib=None
		self._mass=1
		self.MgCl2calib=None
		self.MgGconc=1.1
		self.graphing=False
		self.titrations=None


	def sc_ATP(self, df=None, _from='MgCl2', **kwargs): # adtp_col=None, adtp=None, ATP_calib=None, _mass=1, MgCl2calib=None, MgGconc=None, graphing=False):
		'''
		Calibrate the ATP signal
		first calibrates to MgCl2 titration,

		
		
		Requires df with at least 'Event' and column labeled 'ADP' or 'ATP' 
		'''
		
		if df is not None:pass
		else:df=self.df
		self.df=df

		# A few checks
		for k, v in kwargs.items():
			if k in ['adtp_col', 'adtp', 'ATP_calib', 'MgCl2calib', 'MgGconc', 'graphing', '_mass', 'titrations']:
				if v !=None: setattr(self, k, v)

		# Assign default
		MgGalpha=1
		self.calibration={'parameter':self.adtp,
					'MgG': MgGalpha,
					'sc': [[0],[0]],
					'slope':float('nan'),
					'intercept':float('nan'),
					'r2': float('nan'),
					'predicted':[[0]]
					}

		# Sort adtp col
		if self.adtp_col!=None: adtp_col=self.adtp_col
		else:adtp_col=[c for c in df.columns if 'ADP' in c or 'ATP' in c][0]

		if type(self.MgGconc) is float:
			# Calibrate to MgCl2
			# Needs MgG as event
			# MgG concentration 1.1uM by default
			try:
				MgG_idx=np.where(df.index==df.loc[df['Event']=='MgG'].index[0])[0][0]
				bsln=df[adtp_col].iloc[(MgG_idx-6):(MgG_idx-2)].mean()
				MgG=df[adtp_col].iloc[(MgG_idx+2):(MgG_idx+6)].mean()
				MgGalpha=self.MgGconc/(MgG-bsln)
			except Exception as e:
				print(f'MgG calibration ERROR: {e}')

		if type(self.MgCl2calib) is dict:
			# Calibrate to MgCl2 titration
			# Requires a dictionary {MgCl2: concentration}
			mgsc={} # store mean calibration values to fit sc as fluo: [MgCl2]
			for t, conc in self.MgCl2calib.items():
				fluo=df.loc[df['Event']==t,'ATP']
				fluo=fluo.iloc[int(0.25*len(fluo)):int(0.75*len(fluo))].mean() # Select 50% in middle and average
				mgsc.update({fluo:conc})
			mgsc=pd.DataFrame.from_dict(mgsc, orient='index').sort_index()

			# Calibrate the signal
			slope, intercept, r2, predicted = linear_fit(mgsc.index.values, mgsc.values)
			self.mgcl2_calib = {'sc_raw':mgsc,
							'sc': predicted,
							'slope':slope,
							'intercept':intercept,
							'r2': r2,
							'predicted':predicted}

			self.df['MgCl2']= self.mgcl2_calib['slope']*self.df['Amp_bckg_corrected'] + self.mgcl2_calib['intercept']

		# Differenciate between ATP and ADP calib
		# At this stage, it returns the calibrated signal after
		# the first calibration
		# If no 'ADP' or 'ATP' event in place,
		# returns the non calibrated signal
		
		if True:

			if _from == 'MgCl2': adtp_col=_from


			# only select events that are present in ATP_calib
			self.ATP_calib={k:v for k,v in self.ATP_calib.items() if k in self.titrations}

			if True:
				
				sc={} # as x:y ([MgCl2]: [A(D/T)P])
				for event, conc in self.ATP_calib.items():
					if conc==0:_idx=0.9 	# chose juste before first titration
					else: _idx=0.1 			# otherwise 10-15%% after just after
					fluo=self.df.loc[df['Event']==event, adtp_col]
					fluo=fluo.iloc[int((_idx*len(fluo))):int((_idx+0.05)*len(fluo))].mean()
					sc.update({fluo: conc})
				
				# Retrieve linear sc
				x=[k for k,_ in sc.items()]
				y=[v for _,v in sc.items()]

				# Get standard curve
				slope, intercept, r2, predicted = linear_fit(x, y)

				# This is for ADP calib.
				# Seems to pick up wrong 0
				# so only do 2 points cal
				if r2<=0.9:
					x, y = x[:-1], y[:-1]
					slope, intercept, r2, predicted = linear_fit(x, y)


				sc_raw=[x,y]# store to check sc

				# Keep calibration for future analysis
				self.calibration={
					'parameter': self.adtp,
					'MgG': MgGalpha,
					'mgcl2_calib': self.mgcl2_calib,
					'sc_raw':sc_raw,
					'sc': predicted,
					'slope':slope,
					'intercept':intercept,
					'r2': r2
					}

				if self.graphing:
					plt.scatter(sc_raw[0], sc_raw[1])
					plt.plot(X, predict)
					#plt.ylim(0,6)
					#plt.xlim(0,6)
					plt.show()

				return self.calibration
			else: pass
		# except Exception as e:
		# 	print(f'ATP Calib ERROR:{e}')
		
		# If any error, return the raw ATP column
		return self.calibration


	def calibrate(self, df=None, _method='init', correct_titration=False, with_slope=True, **kwargs):
		'''
		Wrapper function to calibrate the signal based
		By default, this expects MgCl2 calibration and ATP or ADP calibration within the experiment.
		if the _method is altered, expects slope and intercept to calibrate the ATP Amp signal.
		Expects a df that has 'ATP' column at least.
		Any kwarg argument is passed to the sc_ATP() function of this ATP class
		'''

		# A few checks
		if df is not None:pass
		else:df=self.df
		self.df['Amp']=df['ATP']
		for k, v in kwargs.items():
			if k in ['adtp_col', 'adtp', 'ATP_calib', 'MgCl2calib', 'MgGconc', 'graphing', '_mass', 'correct_titration']:
				if v !=None: setattr(self, k, v)

		if correct_titration:
			#Expects a list of [[previous event, event to correct], [previous event, event to correct],...]
			if 'titration' in kwargs: self.titrations=kwargs['titrations']
			self.df['Amp_bckg_corrected']=self.delete_substrate_background(pd.concat([df['ATP'], df['Event']], 
																axis=1),
														correct_titration,
														self.titrations)

		# This calibrates to MgCl2
		sc=self.sc_ATP(df=self.df, adtp_col='Amp_bckg_corrected', _from='MgCl2')

		# Calibrate ATP from the ATP titration
		if _method == 'init':
			self.df['ATP']=sc['slope']*self.df['MgCl2'] + sc['intercept']

		return sc


	def delete_substrate_background(self, df, correct_titration, titrations):
		'''
		Delete the shift in fluorescence due to chemical titrations.
		Levels the 90-95% pre and 5-10% post.
		Requires correct_titration, list of titration that needs correction, and assay titration
			to retrieve the previous titration from the one that needs correction.
		'''

		# First, create the titration pairs
		titration_pairs=[]
		for t in correct_titration:
			for i in range(len(titrations)):
				if titrations[i]==t:
					titration_pairs.append([titrations[i-1],t])


		for t in titration_pairs:
			prev_sig=df.loc[df['Event']==t[0]]['ATP']
			prev_sig=prev_sig.iloc[int(len(prev_sig)*0.9):int(len(prev_sig)*0.95)].mean()
			post_titr=df.loc[df['Event']==t[1]]['ATP']
			_idx=post_titr.index[0]-df.index[0]-2
			post_titr=post_titr.iloc[int(len(post_titr)*0.05):int(len(post_titr)*0.1)].mean()
			df.iloc[_idx:,0]=df['ATP'].iloc[_idx:]-(post_titr-prev_sig)

		return df['ATP']


def hill(S, Vmax, Km, n=1):
	S=np.abs(S)
	return (Vmax*S**n)/(Km+S**n)


def fit_hill_curve(x, y):
	'''
	
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



def sort_events(df):
	"""
	Distributes events name and comments (text) to respective chambers.
	Takes the yet unprocessed df (i.e. exported csv from DatLab)
	"""
	try:
		df['Event Name'] = df['Event Name'].fillna(method='ffill')
		df['Chamber'] = df['Chamber'].fillna(method='ffill')
		df['A_Event'] = np.where((df['Chamber'] == 'Left') | (df['Chamber'] == 'Both'),df['Event Name'], float('nan'))
		df['A_Event_text'] = np.where((df['Chamber'] == 'Left') | (df['Chamber'] == 'Both'),df['Event Text'], float('nan'))
		df['B_Event'] = np.where((df['Chamber'] == 'Right') | (df['Chamber'] == 'Both'),df['Event Name'], float('nan'))
		df['B_Event_text'] = np.where((df['Chamber'] == 'Right') | (df['Chamber'] == 'Both'),df['Event Text'], float('nan'))
		df['A_Event'] = df['A_Event'].fillna(method='ffill')
		df['A_Event_text'] = df['A_Event_text'].fillna(method='ffill')
		df['B_Event'] = df['B_Event'].fillna(method='ffill')
		df['B_Event_text'] = df['B_Event_text'].fillna(method='ffill')
	except Exception as e: print(f"ERROR: sort_events(): {e}")

	return df



def detect_protocol(expe_df):
	'''
	return the tiration protocol from the event.
	Need to run sort_event first
	'''
	titrations={'A': [],
				'B': []}

	for chb, titr in titrations.items():
		try:
			for i in expe_df[f"{chb}_Event"].to_list():
				if (i not in titr) and (i != float('nan')):
					titr.append(i)
		except Exception as e:
			print(f"detect_protocol() ERROR: {e}")

	return titrations



def split_chamber_df(df, chb, fluo=None):

	# 1) Sort columns names
	kPa=False
	for i in [': O2 pres', ': O2 conc']:
		if 'pres' in i: kPa=True
		try: PO2_col=[c for c in df.columns if f"{chb}{i}" in c][0]
		except:pass

	JO2_col=[c for c in df.columns if f"{chb}: O2 flux" in c][0]

	Amp=False
	if any('Amp' in c for c in df.columns):
		Amp_col=[c for c in df.columns if f"{chb}: Amp " in c][0]
		Amp=True

	pX=False
	if any('pX' in c for c in df.columns):
		pX_col=[c for c in df.columns if f"{chb}: pX " in c][0]
		pX=True

	# 2) Sort the fluorescence
	fluo_col=None
	if fluo is None:
		if Amp: fluo, fluo_col = 'Amp', Amp_col
		elif pX: fluo, fluo_col = 'pX', pX_col
		else: pass

	else:
		try:
			fluo=[k for k, v in FLUORESCENCE.items() if fluo in v][0]
		except Exception as e:
			print(f"class Experiment retrieving fluorescence ERROR: {e}")

		if Amp: fluo_col=Amp_col
		elif pX: fluo_col=pX_col
		else: pass


	# 3) Retrieve chamber df:
	loc_cols=[PO2_col, JO2_col, f'{chb}_Event']
	if fluo_col != None: loc_cols=loc_cols+[fluo_col]
	chamber_df=df.loc[:,loc_cols].dropna()
	chamber_df=chamber_df.rename(columns=
				{PO2_col:'PO2',
				JO2_col:'JO2',
				f'{chb}_Event': 'Event',
				fluo_col: fluo})

	return chamber_df



def delete_reox(df, o2col='PO2', parameter='JO2', reox_mode='auto', window=[-5,500], **kwargs):
	'''
	Delete reoxygenation if required
	Please precise '_mode'
	By default, 'auto' will take the [O2] as ref and delete any
		[O2] going up.
		It will also delete portions were JO2 is extreme
		and this may be precised in 'window' variable as [min, max].
	If _mode is 'reox':
		Expects 'reox' Events, which will be deleted until next event.
	'''

	df[o2col]=df[o2col].ewm(span=4).mean()
	if reox_mode == 'auto':
		df=df[df[o2col].diff()<=0]

	elif reox_mode == 'reox':
		df=df[df[f"Event"]!='reox'].sort_index()

	return df.reset_index()



class Chamber():
	'''
	Creates class object for chamber parameter.
	Requires the chamber dataframe extracted from experiment class 
	
	expecting:
							'filename':csv,
							'temperature':temperature,
							'chamber':chb,
							'state':_state,
							'df':cdf,
							'mass':weight,
							'fish':fish,
	'''
	def __init__(self, chamber_dict=False, **kwargs):

		if chamber_dict:
			for k, v in chamber_dict.items():
				setattr(self, k, v)


		else:
			self.df = pd.DataFrame(columns=['PO2', 'JO2', 'Amp', 'pX'], index=[])

			self.PO2 =  self.df['PO2']
			self.JO2 =  self.df['JO2']
			self.Amp = self.df['Amp']
			self.pX = self.df['pX']

			self.filename = str('filename')
			self.chamber = 'A' # will be A or B
			self.temperature = 37 # Default
			self.sample = 'mitochondria'
			self.mass = 0
			self.corhort = None
			self.fluo=None


	def zero_PO2(self, _min=None):
		'''
		calibrate 0%
		'''

		if _min==None: _min=self.df['PO2'].min()
		
		self.df['PO2']=self.df['PO2']-_min



	def calculate_JO2(self, window=20, Kpa_to_pmol=True):
		'''
		Recalculates JO2 from PO2.
		Expects PO2 as series and in kPa.
		Requires weight in mg (float)
		'''

		# Moving average
		self.PO2=self.PO2.rolling(window=int(window/4)).mean()

		# Calculate JO2
		# Expressed in kPa/s/mg
		JO2=((((self.PO2.shift())-self.PO2)*CHAMBER_VOLUME)/(mass*2))*1000

		# Exponential moving average to smooth it all.
		JO2=JO2.ewm(span=window, adjust=False).mean()

		# Convert to pmol/s/mg
		if Kpa_to_pmol is True:
			JO2=JO2*(224.47/20.37)

		self.JO2=JO2
		return JO2




	def average_states(self, _mode='time', **kwargs):
		'''
		Retrieves averages and stores in dataframe
		
		Requires:
		- chamber: df of the chamber with at least JO2 column and time index
			events are meant to be sorted prior to this
		- _mode: 'stable' 
		'''



		# A few checks
		if 'mass' in kwargs: _mass=kwargs['mass']
		else: _mass=self.mass
		
		if 'window' in kwargs: window=kwargs['window']
		else: window=5
		
		if 'parameter' in kwargs:
			if type(kwargs['parameter']) is list: parameter_s=kwargs['parameter']
			else: parameter_s=[kwargs['parameter']]
		else:
			parameter_s= self.df.columns.to_list()
			parameter_s=parameter_s[np.logical_not(np.isnan(my_array))]



		
		favdf=[] # empty list to create final dataframe
		for parameter in parameter_s:

			# Create the row that will store infor ation
			row=pd.DataFrame(index=[parameter], columns=self.titrations)

			# Select parameter only
			self.df[parameter].astype('float64')
			pdf=pd.DataFrame(self.df.loc[:,[parameter,'Event']]).reset_index()

			for event, pedf in pdf.groupby(pdf['Event']):
				_mean=float('nan')
				if event in self.titrations:

					# First sort time in case of reoxygenation deletions
					# and smooth signal
					pedf.loc[:,parameter]=pedf.loc[:,parameter].sort_index().rolling(window=window, min_periods=1).mean()
					
					_idx=int(len(pedf)*0.25)
					if _mode == 'stable':
						# If run into anoxia, only select start of succinate
						select=pedf.iloc[_idx:len(pedf)-int(_idx/2)]
						if ('anoxia' in self.protocol) and (event=='S'):
							select=pedf.iloc[_idx:_idx+30]
						
						if '_exclude' in kwargs:_exclude=kwargs['_exclude']
						else:_exclude=0.1 #Arbitrary value that seems to work
						
						# Select stable part:
						_fder=pedf[parameter].diff()/pedf['index'].diff()
						pedf.loc[:,'fsder']=_fder.values
						# delete parts with high/low fsder
						mask=(pedf['fsder']<_exclude)&(pedf['fsder']>-_exclude)
						select=pedf[parameter].loc[mask]

					else: 	
						# Select only 60% in between
						select=pedf[parameter].iloc[_idx:len(pedf)-_idx]

					# Calculate mean and averages
					# Note that CCCP, maximum is taken into account
					_mean=select.mean()
					_sd=select.std()

					# ============= Parameter specific:

					if parameter =='JO2':
						if event in ['CCCP', 'TMPD-Asc','FCCP']:
							_mean=select.max()
						if event in ['Oli', 'KCN', 'Azide']:
							_mean=select.min()



					# Check if _sd , 10% of _mean
					_error=_sd/_mean
					if _error>=0.5:
						print(f"/!| Large variation in {event} ; mean={_mean} sd={_sd}|!|")
						pass

					row.loc[parameter,event]=_mean
					#row.loc[f"{parameter}_sd",event]=_sd

			favdf.append(row)
		favdf=pd.concat(favdf)
		return favdf



class Experiment():
	'''
	Parameters:
		- csv_path (str): csv path in full
		- Optional:
			- 'fluo': 
			- 'fluo_channel':
	'''

	def __init__(self, csv_path, fluo=None, **params):


		# --- sort parameters
		self.filename=csv_path.split('/')[-1]

		self.__defaultParams__()
		for k, v in params.items():
			try:
				if k in ['sample', 'subsample', 'corhort', 'mass', 'unit']:
					if type(v) is not dict:
						v={'A':v, 'B':v}
				setattr(self, k, v)
			except Exception as e:
				print(f"Experiment class init ERROR: {e}")
		#self.params=self.__dict__

		# --- Sort fluorescence 
		if type(fluo) is not dict:
			fluo={'A': fluo, 'B': fluo}
		self.fluo=fluo


		# --- sort main df
		try:
			self.df_raw=sort_events(pd.read_csv(csv_path,
								low_memory=False,
								encoding= 'unicode_escape'))
		except Exception as e:
			print(f"Experinment class init ERROR: {e}")
			print("/! Cannot create experiment. ABORD.")
			return


		self.titrations=detect_protocol(self.df_raw)

		self.split_chambers()	



	def __defaultParams__(self):

		self.temperature = 20
		self.stirrer_speed=750
		self.POS_V=800
		self.pressure = 1 #atm
		self.recording_interval = 2

		self.sample = {'A': None, 'B': None}
		self.subsample = {'A': None, 'B': None}
		self.corhort = {'A': None, 'B': None}
		self.mass = {'A': 0, 'B': 0}
		self.unit = {'A': 'mg', 'B': 'mg'}



	def split_chambers(self):

		chambers=[]

		for chb in ['A', 'B']:
			
			chamber={'chamber': chb}
			chamber['fluo']=self.fluo[chb]
			chamber['titrations']=self.titrations[chb]
			chamber['df']=split_chamber_df(self.df_raw, chb, fluo=chamber['fluo'])

			# add params of the whole experiment

			for k, v in self.__dict__.items():
				if k not in list(chamber.keys()):
					if type(v) is dict:
						chamber.update({k:v[chb]})
					else:
						chamber.update({k:v})

			chambers.append(chamber)

			#print(chamber)
		# --- Sort parameters for each chambers
		# A, B = self.params, self.params
		# for k, v in self.params.items():
		# 	if k in ['sample', 'subsample', 'corhort', 'mass', 'unit']:
		# 		if type(v) is dict:
		# 			A.update({k:v['A']})
		# 			B.update({k:v['B']})

		# A['fluo'], B['fluo'] = self.fluo['A'], self.fluo['B']
		# A['titrations'], B['titrations'] = self.titrations['A'], self.titrations['B']

		# A['df']=split_chamber_df(self.df_raw, 'A', fluo=A['fluo'])
		# B['df']=split_chamber_df(self.df_raw, 'B', fluo=B['fluo'])

		self.A=Chamber(chambers[0])
		self.B=Chamber(chambers[1])



def background_titrations():
	'''
	Chamber A: ATP
	Chamber B: ADP

	Events as chemical cc. e.g. 'ADP 0.25'
	'''
	ADENYLATES=['ADP', 'ATP']
	STANDARD={'ADP':{}, 'ATP':{}}
	temperatures=[18,20,24,26,27,28,30,32,34]
	for a in ADENYLATES:
		STANDARD[a].update({t:pd.DataFrame() for t in temperatures})


	OUTLIERS=['ATP_ADP_Calibration_Trial_T2_34_16.6.23.csv',
			'ATP_ADP_Calibration_Trial_T4_34_16.6.23.csv']
	
	experiments=[] #Storage after csv extraction
	# I - extract csvs
	subfolders= [f"{DATAPATH}Background Temperature Calibrations/Background {t}/" for t in temperatures]
	for subfolder in subfolders:
		_csvs=[f"{subfolder}{f}" for f in os.listdir(subfolder) if '.csv' in f]
		for csv in _csvs:
			filename=csv.split('/')[-1]
			temperature=int(filename.split('_')[-2])
			trial_nb=filename.split('_')[-3]

			params={"temperature":temperature,
				'filename':filename,
				'trial':trial_nb}

			expe=Experiment(csv, fluo='ATP', **params)
			experiments.append(expe)
			print(f"extracted {filename}")

	# II - analyse each experiment
	for e in experiments:
		for chamber in [e.A, e.B]: #Avoid writing things twice

			print(f"Processing {chamber.filename} - {chamber.chamber}")
			#print(chamber.titrations)

			if chamber.chamber=='A':adenylate='ATP'
			else: adenylate='ADP'

			# 1) calibrate each chamber: MgG=0, MgCl2=125uM
			mgsc={} # store sc for magnesium
			for t in ['MgG', 'MgCl2']:
				_t=chamber.df.loc[chamber.df['Event']==t,'ATP']
				_t=_t.iloc[int(0.25*len(_t)):int(0.75*len(_t))].mean() # Select 50% in middle and average
				mgsc[t]=_t
			# Calibrate the signal
			chamber.df['ATP']= (chamber.df['ATP']-mgsc['MgG'])*(125/(mgsc['MgCl2']-mgsc['MgG']))

			# 2) Extract titration means
			sc={0:mgsc['MgCl2']} # store the whole standard curve
			for titr in chamber.titrations[1:]:
				if 'A' in titr: # dicard MgG and MgCl2
					conc=float(titr.split(' ')[-1])
					_t=chamber.df.loc[chamber.df['Event']==titr,'ATP']
					val=_t.iloc[int(0.25*len(_t)):int(0.75*len(_t))].mean() # Select 50% in middle and average
					sc[conc]=val

			sc=pd.DataFrame.from_dict(sc, orient='index', columns=[chamber.trial])
			sc=sc-sc.iloc[0]
			chamber.sc=sc

			# Append to dictionary
			if ((chamber.filename in OUTLIERS) and (adenylate=='ADP')):pass

			else:STANDARD[adenylate][chamber.temperature]=pd.concat([STANDARD[adenylate][chamber.temperature], sc], ignore_index=False, axis=1)
	
	# # Graphing
	# fig, ax = plt.subplots(len(ADENYLATES), len(temperatures))
	# for s in range(len(ADENYLATES)):
	# 	for t in range(len(temperatures)):
	# 		for c in list(STANDARD[ADENYLATES[s]][temperatures[t]].columns):
	# 			ax[s,t].scatter(STANDARD[ADENYLATES[s]][temperatures[t]].index, STANDARD[ADENYLATES[s]][temperatures[t]][c].values)
	# 			ax[s,t].set_title(f"{ADENYLATES[s]} - {temperatures[t]}ºC")
	# 			ax[s,t].set_ylim(-2000,100)
	# 			ax[s,t].set_yticklabels([])
	

	# mng = plt.get_current_fig_manager()
	# mng.full_screen_toggle()
	# plt.subplots_adjust(wspace=0)
	# plt.show()


	# Extract linear fit for titration. MgCl2 - A(D/T)P
	fit=STANDARD
	for a in ADENYLATES:
		for t in temperatures:

			STANDARD[a][t]=STANDARD[a][t].fillna(STANDARD[a][t].rolling(window=3, min_periods=1).mean())

			#Get standard curve
			lm=linear_model.LinearRegression()
			x=np.array(STANDARD[a][t].index).reshape(-1,1)
			y=STANDARD[a][t].mean(axis=1).values

			model=lm.fit(x,y)
			r2=lm.score(x,y)
			slope=lm.coef_[0]
			intercept=lm.intercept_

			X=np.arange(0, 5.1, 0.1)
			Y=model.predict(X.reshape(-1,1))
			
			fit[a][t]=pd.DataFrame(Y, index=X)


	# Create the ADP/ATP table with constant adenylate
	adp_atp=pd.DataFrame(index=X)
	adp_atp_ratio=pd.DataFrame(index=X)

	fig, ax=plt.subplots(2,len(temperatures))
	for t in range(len(temperatures)):
		adp=fit['ADP'][temperatures[t]]
		atp=fit['ATP'][temperatures[t]].iloc[::-1].set_index(adp.index)
		# Make it all positive
		_min=min(list(adp.min()), list(atp.min()))[0]
		adp=adp-_min
		atp=atp-_min

		_ratio=atp/adp

		adp_atp[f"{temperatures[t]}_ADP"]=adp
		adp_atp[f"{temperatures[t]}_ATP"]=atp
		adp_atp_ratio[temperatures[t]]=_ratio

		# do the polyfit
		x=_ratio.index.values
		poly=np.polyfit(x, _ratio[0].values, 2)
		predicted=(poly[0]*(x**2))+(poly[1]*x)+poly[2]
		#Vmax, pP50, _hill, predicted, r2=fit_hill_curve(_ratio[0].values, _ratio.index.values)

		ax[0,t].scatter(adp.index, adp, c='g')
		ax[0,t].scatter(atp.index, atp, c='r')
		ax[0,t].set_ylim(0,1500)
		ax[0,t].set_title(temperatures[t])
		ax[1,t].scatter(_ratio.index, _ratio)
		ax[1,t].plot(predicted)

		ax[1,t].set_xlim(0,5)

	print(adp_atp_ratio)
	#adp_atp_ratio.to_csv("ratio.csv")
	adp_atp.to_csv("adp_atp.csv")
	plt.show()



def exctract_PW_csv(data_path=f"{current_path}/ATP_ADP_PW_2022/", repo_path=REPO_PATH, repo_tabname='Parrot_Wrasse_Chamber_Tisue', **kwargs):
	
	experiments=[] # Main store

	# 1) Extract mass and other relevant information from the log
	repo=pd.read_excel(repo_path, sheet_name=repo_tabname)
	repo=repo.loc[:,['ATP_ADP_Experiments','Chamber_A','Chamber_B']]

	# 2) Exctract csv 
	for temperature in kwargs['temperatures']:
		_csvs=[c for c in os.listdir(f"{data_path}{temperature}/Exported/") if 'csv' in c]
		for csv in _csvs:
			try:
				print(f"Exctracting {csv}...")
				csv_path=f"{data_path}{temperature}/Exported/{csv}"

				# repo does not have 'Exported' in the name so add it to retrieve params:
				if 'Exported' in csv:
					_rename=csv[:-4].split('_')
					_rename=f"{_rename[0]}_{_rename[1]}_{_rename[2]}_{_rename[4]}_{_rename[5]}_{_rename[6]}"

				# Retrieve parameters from repo and arrange to dict
				params=repo.loc[repo['ATP_ADP_Experiments']==_rename]

				params={k:v[0] for k,v in params.to_dict(orient='list').items()}

				params.update({
								'species':'Parrot',
								'protocol':'anoxia_steady',
								'temperature':temperature,
								'mass':{'A':params['Chamber_A'],
										'B':params['Chamber_B']
										},
								'corhort':{'A':'ATP',
										'B':'ADP'
										}
								})

				# Create experiment object
				expe=Experiment(csv_path, fluo='ATP', **params)
				experiments.append(expe)
				plt.plot(expe.A.df['JO2'])
				plt.plot(expe.A.df['ATP'])
				plt.show()
				plt.clear()

				plt.plot(expe.B.df['JO2'])
				plt.plot(expe.B.df['ATP'])
				plt.show()
			except Exception as e: print(f"ERROR: exctracting_csv(): {csv}: {e}")
	return experiments



def exctract_BW_Steady_csv(data_path=f"{current_path}/Banded_Wrasse_ATP_SteadyState/", repo_path=REPO_PATH, repo_tabname='ATP Steadystate Tissue Weights', **kwargs):
	
	experiments=[] # Main store

	# 1) Extract mass and other relevant information from the log
	repo=pd.read_excel(repo_path, sheet_name=repo_tabname)
	repo=repo.loc[:,['File Name','Chamber A','Chamber B']]

	# 2) Exctract csv 
	for temperature in kwargs['temperatures']:
		_csvs=[c for c in os.listdir(f"{data_path}{temperature}/Exported/") if 'csv' in c]
		for csv in _csvs:
			try:
				print(f"Exctracting {csv}...")
				csv_path=f"{data_path}{temperature}/Exported/{csv}"

				_rename=csv[:-4]

				# Retrieve parameters from repo and arrange to dict
				params=repo.loc[repo['File Name']==_rename]

				params={k:v[0] for k,v in params.to_dict(orient='list').items()}

				params.update({
								'species':'Banded',
								'protocol':'steady',
								'temperature':temperature,
								'mass':{'A':params['Chamber A'],
										'B':params['Chamber B']
										},
								'corhort':{'A':'ATP',
										'B':'ADP'
										}
								})

				# Create experiment object
				expe=Experiment(csv_path, fluo='ATP', **params)
				experiments.append(expe)
			except Exception as e: print(f"ERROR: exctracting_csv(): {csv}: {e}")
	return experiments



def exctract_BW_anoxia_csv(data_path=f"{current_path}Banded_Wrasse_ATP_Anoxia/", repo_path=REPO_PATH, repo_tabname='ATP Rundown Tissue Weights', **kwargs):
	
	experiments=[] # Main store

	# 1) Extract mass and other relevant information from the log
	repo=pd.read_excel(repo_path, sheet_name=repo_tabname)
	repo=repo.loc[:,['File Name','Chamber A','Chamber B']]

	# 2) Exctract csv 
	for temperature in kwargs['temperatures']:
		_csvs=[c for c in os.listdir(f"{data_path}{temperature}/Exported/") if 'csv' in c]
		for csv in _csvs:
			try:
				print(f"Exctracting {csv}...")
				csv_path=f"{data_path}{temperature}/Exported/{csv}"

				_rename=csv[:-4]

				# Retrieve parameters from repo and arrange to dict
				params=repo.loc[repo['File Name']==_rename]

				params={k:v[0] for k,v in params.to_dict(orient='list').items()}

				params.update({
								'species':'Banded',
								'protocol':'steady',
								'temperature':temperature,
								'mass':{'A':params['Chamber A'],
										'B':params['Chamber B']
										},
								'corhort':{'A':'ATP',
										'B':'ADP'
										}
								})

				# Create experiment object
				expe=Experiment(csv_path, fluo='ATP', **params)
				experiments.append(expe)
			except Exception as e: print(f"ERROR: exctracting_csv(): {csv}: {e}")
	return experiments



def fit_anoxia_curve(experiment, graphing=False):

	print(f"Processing {experiment.filename}...")
	# ========== Fit and exctract profiles
	for chamber in [experiment.A, experiment.B]:
		try:

			# Calibrate PO to 0
			chamber.df['PO2']=chamber.df['PO2']-chamber.df['PO2'].min()
			print(chamber.titrations)


			# Select Anoxia portion
			tdf=chamber.df.loc[(chamber.df['Event']=='S')]

			# Cleanup the portion
			tdf=tdf.iloc[40:,:].loc[:,['PO2','JO2']].ewm(span=10).mean()

			# prepare for curve fitting
			tdf=tdf.set_index('PO2').sort_index()

			# Fit curve to experimental PO2 
			fJO2max, fP50, _hill, fitted, r2 = fit_hill_curve(tdf.index.values, tdf['JO2'].values)

			# Extend curve to 0 -> 100% PO2
			# Standardise PO2 for comaprison between samples
			X=np.arange(0, MAX_kPA, 0.01)
			predicted=hill(X, fJO2max, fP50, _hill)
			predicted=pd.DataFrame(predicted, index=X, columns=['JO2'])

			# Calculate area under the curve from Standardised
			for i in [21, 10, 5, 2.5, 1.25, 0.5]:
				_trimmed_df=predicted.loc[predicted.index<i]
				auc=np.trapz(_trimmed_df['JO2'], _trimmed_df.index, 0.1)

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

			if graphing:
				plt.scatter(tdf.index.values, tdf['JO2'].values, c='r')
				plt.scatter(X, predicted, c='b')
				plt.show()

		except Exception as e:
			print(f"ERROR: fitting_curve_main - {chamber.filename}: {e}")



def main(save_sc=False):

	experiments=[]
	standard_curves=[]
	averages=[]
	bck_sc={'ADP':{}, 'ATP':{}}

	# Open Background titration curves
	bckg_df=pd.read_csv(f"{current_path}/Jules_Analysis_CSVs/adp_atp.csv", index_col=0)
	for c in bckg_df.columns:
		temp=int(c.split('_')[0])
		aden=c.split('_')[1]
		slope, intercept, r2, predicted = linear_fit(bckg_df.index, bckg_df[c].values)
		if slope>0: slope=-slope
		bck_sc[aden].update({temp: [slope, intercept]})
	for a in ['ADP', 'ATP']:
		bck_sc[a][25]=[(bck_sc[a][24][0]+bck_sc[a][26][0])/2, (bck_sc[a][24][1]+bck_sc[a][26][1])/2]


	# Open ratio df
	ratios=pd.read_csv(f"{current_path}/Jules_Analysis_CSVs/ADTP_total.csv", index_col=0)

	# Exctract csvs
	PW_expe=exctract_PW_csv(parameters={'fluo':'ATP'}, temperatures=[24,26,28,30,32,34])
	BW_steady_expe=exctract_BW_Steady_csv(parameters={'fluo':'ATP'}, temperatures=[18,20,24,25,26,27,30])
	BW_anoxia_expe=exctract_BW_anoxia_csv(parameters={'fluo':'ATP'}, temperatures=[18,20,26,27,28,30])
	# 
	experiments+=PW_expe
	experiments+=BW_steady_expe
	experiments+=BW_anoxia_expe


	# Process each assay 
	for experiment in experiments:
		if experiment.filename not in OUTLIERS:

			#try:

			# Define calibrations and other things depending on the species specific experiments:
			if experiment.species == 'Parrot':
				MgCl2_calib={'CaCl2': 0, #mM
							'MgCl2_1':62.5,
							'MgCl2_2':125}

				correct_titration=['PMG','S']

			elif experiment.species == 'Banded':
				MgCl2_calib={
							'Ca2+': 0, #mM
							'MgCl2_1':62.5,
							'MgCl2_2':125
							}
				if experiment.protocol=='steady':
					correct_titration=['PMG','S','Oli']
				elif experiment.protocol=='anoxia':
					correct_titration=['PMG','S']
				

			# Loop through each chamber
			print(f"Processing {experiment.filename}...")
			for chamber in [experiment.A, experiment.B]:

				# ====================== ATP CALIBRATION ===============================
				# Calibrate ATP and ADP from the chamber calibration itself

				# First, calibrate ATP signal and retrieve calibration data
				# Preload parameters
				if chamber.chamber=='A':
					adtp='ATP'
					adtp_calib=ATP_calib

				else:
					adtp='ADP'
					adtp_calib=ADP_calib
					correct_titration+=['Bleb', 'Ouab']

				# Prepare parameters for ATP calibration
				atp_params={
							'df':chamber.df,
							'_mass':float(chamber.mass),
							'adtp':adtp,
							'adtp_col':'ATP',
							'ATP_calib':adtp_calib,
							'MgCl2calib':MgCl2_calib,
							'graphing':False,
							'titrations':chamber.titrations
							}
			
				# See Fluorescence.ATP.sc_ATP() for more info
				atp=ATP(**atp_params)
				sc_ATP=atp.calibrate(correct_titration=correct_titration)
				chamber.df=atp.df
				chamber.df['ATP_1']=chamber.df['ATP']
				# Append dict to the standard curve list for latter processing
				standard_curves.append({
					'filename':chamber.filename,
					'temperature':chamber.temperature,
					'adtp':adtp,
					'species':chamber.species,
					'mass':chamber.mass,
					'MgCl2_slope':sc_ATP['mgcl2_calib']['slope'],
					'MgCl2_intercept':sc_ATP['mgcl2_calib']['intercept'],
					'MgCl2_r2':sc_ATP['mgcl2_calib']['r2'],
					f'{adtp}_slope':sc_ATP['slope'],
					f'{adtp}_intercept':sc_ATP['intercept'],
					f'{adtp}_r2':sc_ATP['r2'],
					})
				# ----------------------------- End of Standard curve -------------------------------


				# =========================== Calibrate on bckgd titrations =========================
				slope, intercept = bck_sc[adtp][chamber.temperature][0], bck_sc[adtp][chamber.temperature][1]
				# [ATP] corrected from [MgCl2]
				chamber.df['ATP_2'] = (chamber.df['MgCl2']-intercept)/slope
				# ---------------------------- End of bckgd calibration ------------------------------
				

				# =========================== Calibrate on MgCl2 and ratio =========================
				# Retrieve adtp total profile from the ratio linear regression:
				if adtp=='ADP':
					Y=ratios[str(chamber.temperature)].values
				elif adtp=='ATP':
					# Have to reverse the difference
					Y=-ratios[str(chamber.temperature)][::-1].values
				
				# corrects for MgCl2 experimental variation / re-calibrate from ATP calibration
				# This is the difference between [MgCl2] from titration curves to the experimental [MgCl2] on ADP or ATP addition
				b=max(Y)-sc_ATP['mgcl2_calib']['intercept']

				slope, intercept, r2, predicted = linear_fit(Y, ratios.index)
				chamber.df['ATP_3']=slope*(chamber.df['MgCl2']+b)+intercept
				# ---------------------------- End of bckgd calibration ------------------------------


				# Get the rates
				# rate = difference * 2mL / (2sec * mass)
				for c in [1,2,3]:
					chamber.df[f'ATP_flux_{c}']=((chamber.df[f'ATP_{c}'].ewm(span=5).mean().diff(periods=5))/chamber.mass).ewm(span=5).mean()*1000
					chamber.df[f'PO_{c}']=(chamber.df[f'ATP_flux_{c}']*2/chamber.df['JO2']).ewm(span=5).mean()
					

				# =============== Set boundaries ============
				chamber.df.loc[(chamber.df['JO2']<=0)|(chamber.df['JO2']>=300)]=np.nan
				for c in [1,2,3]:
					if chamber.species=='Parrot':
						chamber.df[f'ATP_{c}']=chamber.df[f'ATP_{c}']*100
						if chamber.temperature <=32:
							chamber.df.loc[:,f'ATP_flux_{c}']=abs(chamber.df.loc[:,f'ATP_flux_{c}'])
					
					else:
						chamber.df[f'ATP_{c}']=chamber.df[f'ATP_{c}']*1000
						if chamber.temperature <=27:
							chamber.df.loc[:,f'ATP_flux_{c}']=abs(chamber.df.loc[:,f'ATP_flux_{c}'])

					chamber.df.loc[(chamber.df[f'ATP_{c}']<=0)|(chamber.df[f'ATP_{c}']>=100),[f'ATP_{c}']]=np.nan
					chamber.df.loc[(chamber.df[f'PO_{c}']<=-1)|(chamber.df[f'PO_{c}']>=5), [f'PO_{c}']]=np.nan

				
				# ================================ Average states ===================================
				parameters=['JO2']
				parameters+=[f'ATP_flux_{c}' for c in [1,2,3]]
				parameters+=[f'PO_{c}' for c in [1,2,3]]
	
				favdf=chamber.average_states(parameter=parameters, _mode='time')
				for k,v in chamber.__dict__.items():
					if k in ['filename', 'temperature', 'mass', 'species', 'protocol', 'chamber','fish']:
						favdf[k]=[v]*len(favdf)
				averages.append(favdf)


				# fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True)
				# x=chamber.df.index
				# ax1.scatter(x, chamber.df['JO2'].values, c='red')
				# ax2.scatter(x, chamber.df['ATP_3'].values, c='blue')
				# ax3.scatter(x, chamber.df['ATP_rate_3'].values, c='green')
				# plt.show()

	# Concatenate all averages into one df
	favdf=pd.concat(averages, axis=0).dropna(how='all', axis=1)

	# ================= Calculate ratios etc ===========
	favdf.loc[favdf['chamber']=='A', 'CI_pct'] = favdf.loc[favdf['chamber']=='A','ATP_2']/favdf.loc[favdf['chamber']=='A','S']
	favdf.loc[favdf['chamber']=='B', 'CI_pct'] = favdf.loc[favdf['chamber']=='B','ADP_2']/favdf.loc[favdf['chamber']=='B','S']
	favdf['L/P']=favdf['PMG']/favdf['S']


	# ================ Save to csv ======================
	favdf.to_csv("Stable_state_average.csv")
	print(favdf)

	if save_sc==True:
		standard_curves=pd.DataFrame(standard_curves)
		standard_curves.to_csv('Standard_Curves.csv')


def to_prism():
	'''
	Transpose necessary values for Prism
	Require csv of averages rows being parameters and columns criteria and states
	'''
	temp=[]
	vals=[]

	state='CI_pct'

	df=pd.read_csv("Stable_state_average.csv")

	# # # ================ For the states in rows, columns temperature sub col repeats
	# _df=df.loc[(df.parameter=='JO2')&(df.species=='Banded')&(df.chamber=='B'),['temperature','filename']+states].sort_values(['temperature','filename'])
	# print(_df.T)
	# _df.T.to_csv('temp_forPrism.csv')

	# ================ For template Parrot + Banded, temperature rows, cols wrass species and sub col wrasse nb
	# Select parameter first
	_df=df.loc[(df.parameter=='JO2')&(df.species=='Banded')&(df.chamber=='B'),['temperature','filename',state]].sort_values(['temperature','filename'])

	# Select and transpose only the state that matters
	for temperature, tdf in _df.groupby(_df['temperature']):
		temp.append(temperature)
		vals.append(tdf[state].T.values)

	tdf=pd.DataFrame(vals, index=temp)
	tdf.to_csv('temp_forPrism.csv')
	print(tdf)




if __name__ == '__main__':
	main()
	
