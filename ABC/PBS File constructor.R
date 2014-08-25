setwd("/scratch/jc227089/PBSfiles/")

tt<-1:200

script.file<-'/home/jc227089/Quoll_model/ABC/Quoll_model_ABC.R'

for (ii in 1:length(tt)) {
	##create the sh file
	zz = file(paste('jobrun_',tt[ii],'.sh',sep=''),'w')
		cat('##################################\n',file=zz)
		cat('cd $PBS_O_WORKDIR\n',file=zz)
		cat('source /etc/profile.d/modules.sh\n',file=zz) ### Runs an .sh file which allows modules to be loaded
		cat('module load R\n', file=zz)
		cat("R CMD BATCH --no-save --no-restore '--args sr=5000 file.ID=",tt[ii],"' ",script.file,' jobrun_',tt[ii],'.Rout \n',sep='',file=zz)
		cat('##################################\n',file=zz)
	close(zz)
			
	#submit the job
	system(paste('qsub -l nodes=1:ppn=1 -l pmem=500mb -l walltime=03:30:00 jobrun_',tt[ii],'.sh',sep=''))
}
