###build metadata for cag regressions from teddy data
##20191918
#Tierney

birth_data={}
with open('birth_variables.csv') as f:
	for line in f:
		line=line.rstrip().split(',')
		birth_data[line[0]]=','.join(line[1:])


diabetes_data={}
with open('diabetes_variables.csv') as f:
	for line in f:
		line=line.rstrip().split(',')
		diabetes_data[line[0]]=','.join(line[1:])

diabetes_data2={}
with open('t1d_metadata_from_old_analysis.csv') as f:
	for line in f:
		line=line.rstrip().split(',')
		diabetes_data2[line[0]]=','.join(line[1:])

#for each line in sample dataframe
output=[]
with open('sample_level_data.csv') as f:
	for i,line in enumerate(f):
		if i==0:
			continue
		output_line=[]
		line=line.rstrip().split(',')
		if line[5]=='':
			continue
		sample=line[5]
		subject=line[0]
		try:
			birth_info=birth_data[line[0]].split(',')
			diabetes_info=diabetes_data[line[0]].split(',')
			diabetes_info2=diabetes_data2[line[0]].split(',')
		except:
			continue
		age=line[3]
		output_line.append(sample)
		output_line.extend(line[:2])
		output_line.append(age)
		output_line.extend(birth_info[:2])
		output_line.extend(birth_info[4:])
		diabetes_age=diabetes_info[1]
		seroconversion_age=diabetes_info[0]
		if seroconversion_age !='':
			if int(seroconversion_age)<int(age):
				output_line.append(0)
			if int(seroconversion_age)>=int(age):
				output_line.append(1)
			output_line.append(1)
		else:
			output_line.extend([0,0])
		if diabetes_age !='':
			if int(diabetes_age)<int(age):
				output_line.append(0)
			if int(diabetes_age)>=int(age):
				output_line.append(1)
		else:
			output_line.append(0)
		output_line.extend(diabetes_info[2:])
		output.append('\t'.join([str(x) for x in output_line])+'\n')


output=['\t'.join(['sampleID','subjectID','country', 'age', 'mom_age_at_birth', 'gestational_diabetes', 'gender', 'HLA_risk_class', 'location', 'seroconverted_at_sampling','seroconverted_ever','diabetes_at_sampling', 't1d_ever','IAA', 'GADA', 'IA2A', 'ZNT8A', 'ICA'])+'\n']+output
#make sure wgs samples
with open('diabimmune_metadata.tsv','w') as w:
	for line in output:
		w.write(line.replace('FALSE','0'))
