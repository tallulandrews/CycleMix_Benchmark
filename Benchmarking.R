source("/home/tandrew6/scripts/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"


Calculate_accuracy_trio <- function(pred, truth) {
	tab <- table(pred, truth)
	if (nrow(tab) < 3) {
		for (i in c("G1", "S", "G2M")) {
			if (!i %in% rownames(tab)) {
				t_ns <- rownames(tab)
				tab <- rbind(tab, rep(0, ncol(tab)));
				rownames(tab) <- c(t_ns, i)
			}
		}
	}
	if ("None" %in% rownames(tab)) {
		if ("G1" %in% rownames(tab)) {
			tab["G1",] <- tab["G1",] + tab["None",]
		} else {
			tmp <- rownames(tab)
			tmp[tmp=="None"] <- "G1"
			rownames(tab) <- tmp
		}
		tab <- tab[rownames(tab) != "None",]
	}
	if ("G1S" %in% rownames(tab)) {
		tmp <- rownames(tab)
		tmp[tmp=="G1S"] <- "S"
		rownames(tab) <- tmp
	}
	if ("G2" %in% colnames(tab)) {
		tmp <- colnames(tab)
		tmp[tmp=="G2"] <- "G2M"
		colnames(tab) <- tmp
	}
	if (nrow(tab) < 3) {
		for (i in c("G1", "S", "G2M")) {
			if (!i %in% rownames(tab)) {
				t_ns <- rownames(tab)
				tab <- rbind(tab, rep(0, ncol(tab)));
				rownames(tab) <- c(t_ns, i)
			}
		}
	}
	tab <- tab[,c("G1", "G2M", "S")]
	TPs <- c(tab["G1", "G1"], tab["G2M", "G2M"], tab["S","S"])
	FPs <- c(tab["G1", "G2M"]+tab["G1","S"], 
		 tab["G2M", "G1"]+tab["G2M","S"], 
		 tab["S","G1"]+tab["S","G2M"])
	FNs <- c(tab["G2M", "G1"]+tab["S","G1"],
                 tab["G1", "G2M"]+tab["S","G2M"],
                 tab["G1","S"]+tab["G2M","S"])
	names(TPs) <- names(FPs) <- names(FNs) <- c("G1", "G2M", "S")
	total <- sum(tab)
	return(list(TPs=TPs, FPs=FPs, FNs=FNs, total=total))
	#return(sum(TPs)/total)
}

# https://en.wikipedia.org/wiki/Phi_coefficient
calc_MCC <- function(accuracy_stats) {
	G1_TN <- (accuracy_stats$total-accuracy_stats$TPs[1]-accuracy_stats$FPs[1]-accuracy_stats$FNs[1])
	G2M_TN <- (accuracy_stats$total-accuracy_stats$TPs[2]-accuracy_stats$FPs[2]-accuracy_stats$FNs[2])
	S_TN <- (accuracy_stats$total-accuracy_stats$TPs[3]-accuracy_stats$FPs[3]-accuracy_stats$FNs[3])
	MCC_G1 <- (accuracy_stats$TPs[1]*G1_TN - accuracy_stats$FPs[1]*accuracy_stats$FNs[1]) / sqrt((accuracy_stats$TPs[1]+accuracy_stats$FPs[1])*(accuracy_stats$TPs[1]+accuracy_stats$FNs[1])*(G1_TN+accuracy_stats$FPs[1])*(G1_TN+accuracy_stats$FNs[1]))
	MCC_G2M <- (accuracy_stats$TPs[2]*G2M_TN - accuracy_stats$FPs[2]*accuracy_stats$FNs[2])/sqrt((accuracy_stats$TPs[2]+accuracy_stats$FPs[2])*(accuracy_stats$TPs[2]+accuracy_stats$FNs[2])*(G2M_TN+accuracy_stats$FPs[2])*(G2M_TN+accuracy_stats$FNs[2]))
	MCC_S <- (accuracy_stats$TPs[3]*S_TN - accuracy_stats$FPs[3]*accuracy_stats$FNs[3])/sqrt((accuracy_stats$TPs[3]+accuracy_stats$FPs[3])*(accuracy_stats$TPs[3]+accuracy_stats$FNs[3])*(S_TN+accuracy_stats$FPs[3])*(S_TN+accuracy_stats$FNs[3]))
	return(c(MCC_G1, MCC_G2M, MCC_S))
}
set.seed(3891)

output_overall <- c();
output_all <- list();
time <- c();

print("Sasagawa")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Sasagawa_with_GT.R")
Sasagawa_meta <- sub("\\/", "", Sasagawa_meta)
Sask_cyclone <- run_cyclone(Sasagawa_fpkms, species="mmus")
Sask_cyclemix <- run_CycleMix(Sasagawa_fpkms, MGeneSets$Cyclone)
Sask_seurat <- run_Seurat(Sasagawa_fpkms, MGeneSets$Cyclone)
acc_cyclone <- Calculate_accuracy_trio(Sask_cyclone$phase, Sasagawa_meta)
acc_seurat <- Calculate_accuracy_trio(Sask_seurat$phase, Sasagawa_meta)
acc_cyclemix <- Calculate_accuracy_trio(Sask_cyclemix$phase, Sasagawa_meta)
acc <- c(sum(acc_cyclone$TPs)/acc_cyclone$total, sum(acc_cyclemix$TPs)/acc_cyclemix$total, sum(acc_seurat$TPs)/acc_seurat$total, n=sum(Sasagawa_meta != "Unknown"))
output_overall <- rbind(output_overall, acc)
time <- rbind(time, c(Sask_cyclone$time, Sask_cyclemix$time, Sask_seurat$time))
output_all[["Sasagawa"]] <- list(cyclone=acc_cyclone, cyclemix=acc_cyclemix, seurat=acc_seurat);

print("Beuttner")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Buttner_with_GT.R")
Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)
acc_cyclone <- Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta)
acc_seurat <- Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta)
acc_cyclemix <- Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta)
acc <- c(sum(acc_cyclone$TPs)/acc_cyclone$total, sum(acc_cyclemix$TPs)/acc_cyclemix$total, sum(acc_seurat$TPs)/acc_seurat$total, n=sum(Buttner_meta != "Unknown"))
output_overall <- rbind(output_overall, acc)
time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))
output_all[["Beuttner"]] <- list(cyclone=acc_cyclone, cyclemix=acc_cyclemix, seurat=acc_seurat);

print("Chen")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Chen_with_GT.R")
Chen_cyclone <- run_cyclone(Chen_counts, species="mmus")
Chen_cyclemix <- run_CycleMix(Chen_counts, MGeneSets$Cyclone)
Chen_seurat <- run_Seurat(Chen_counts, MGeneSets$Cyclone)
acc_cyclone <- Calculate_accuracy_trio(Chen_cyclone$phase, Chen_meta)
acc_seurat <- Calculate_accuracy_trio(Chen_seurat$phase, Chen_meta)
acc_cyclemix <- Calculate_accuracy_trio(Chen_cyclemix$phase, Chen_meta)
acc <- c(sum(acc_cyclone$TPs)/acc_cyclone$total, sum(acc_cyclemix$TPs)/acc_cyclemix$total, sum(acc_seurat$TPs)/acc_seurat$total, n=sum(Chen_meta != "Unknown"))
output_overall <- rbind(output_overall, acc)
time <- rbind(time, c(Chen_cyclone$time, Chen_cyclemix$time, Chen_seurat$time))
output_all[["Chen"]] <- list(cyclone=acc_cyclone, cyclemix=acc_cyclemix, seurat=acc_seurat);

print("Leng")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Leng_with_GT.R")
Leng_cyclone <- run_cyclone(Leng_counts, species="hsap")
Leng_cyclemix <- run_CycleMix(Leng_counts, HGeneSets$Seurat)
Leng_seurat <- run_Seurat(Leng_counts, HGeneSets$Seurat)
acc_cyclone <- Calculate_accuracy_trio(Leng_cyclone$phase, Leng_meta)
acc_seurat <- Calculate_accuracy_trio(Leng_seurat$phase, Leng_meta)
acc_cyclemix <- Calculate_accuracy_trio(Leng_cyclemix$phase, Leng_meta)
acc <- c(sum(acc_cyclone$TPs)/acc_cyclone$total, sum(acc_cyclemix$TPs)/acc_cyclemix$total, sum(acc_seurat$TPs)/acc_seurat$total, n=sum(Leng_meta != "Unknown"))
output_overall <- rbind(output_overall, acc)
time <- rbind(time, c(Leng_cyclone$time, Leng_cyclemix$time, Leng_seurat$time))
output_all[["Leng"]] <- list(cyclone=acc_cyclone, cyclemix=acc_cyclemix, seurat=acc_seurat);

print("MacDavid")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_MacDavid_with_GT.R")
MacDavid_meta[MacDavid_meta == "g0/g1"] <- "G1"
MacDavid_meta[MacDavid_meta == "g2/m"] <- "G2M"
MacDavid_meta[MacDavid_meta == "s"] <- "S"
MacDavid_cyclone <- run_cyclone(MacDavid_counts, species="hsap")
MacDavid_cyclemix <- run_CycleMix(MacDavid_counts, HGeneSets$Seurat)
MacDavid_seurat <- run_Seurat(MacDavid_counts, HGeneSets$Seurat)
acc_cyclone <- Calculate_accuracy_trio(MacDavid_cyclone$phase, MacDavid_meta)
acc_seurat <- Calculate_accuracy_trio(MacDavid_seurat$phase, MacDavid_meta)
acc_cyclemix <- Calculate_accuracy_trio(MacDavid_cyclemix$phase, MacDavid_meta)
acc <- c(sum(acc_cyclone$TPs)/acc_cyclone$total, sum(acc_cyclemix$TPs)/acc_cyclemix$total, sum(acc_seurat$TPs)/acc_seurat$total, n=sum(MacDavid_meta != "Unknown"))
output_overall <- rbind(output_overall, acc)
time <- rbind(time, c(MacDavid_cyclone$time, MacDavid_cyclemix$time, MacDavid_seurat$time))
output_all[["MacDavid"]] <- list(cyclone=acc_cyclone, cyclemix=acc_cyclemix, seurat=acc_seurat);


# All of them are basically at random chance #
print("Madhessian")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Madhessian_with_GT.R")
Madhessian_meta[Madhessian_meta == "G0"] <- "G1"
Madhessian_cyclone <- run_cyclone(Madhessian_counts, species="hsap")
Madhessian_cyclemix <- run_CycleMix(Madhessian_counts, HGeneSets$Seurat)
Madhessian_seurat <- run_Seurat(Madhessian_counts, HGeneSets$Seurat)
acc_cyclone <- Calculate_accuracy_trio(Madhessian_cyclone$phase, Madhessian_meta)
acc_seurat <- Calculate_accuracy_trio(Madhessian_seurat$phase, Madhessian_meta)
acc_cyclemix <- Calculate_accuracy_trio(Madhessian_cyclemix$phase, Madhessian_meta)
#acc <- c(sum(acc_cyclone$TPs)/acc_cyclone$total, sum(acc_cyclemix$TPs)/acc_cyclemix$total, sum(acc_seurat$TPs)/acc_seurat$total, n=sum(Madhessian_meta != "Unknown"))
#output_overall <- rbind(output_overall, acc)
#time <- rbind(time, c(Madhessian_cyclone$time, Madhessian_cyclemix$time, Madhessian_seurat$time))
#output_all[["Madhessian"]] <- list(cyclone=acc_cyclone, cyclemix=acc_cyclemix, seurat=acc_seurat);
###


## Plots! ##
source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")

print("Plotting")
colnames(output_overall) <- c("Cyclone", "CycleMix", "Seurat", "n")
rownames(output_overall) <- c("Sasagawa", "Buettner", "Chen", "Leng", "MacDavid")
colnames(time) <- c("Cyclone", "CycleMix", "Seurat")
rownames(time) <-  c("Sasagawa", "Buettner", "Chen", "Leng", "MacDavid")
output_overall_strerr <- sqrt(output_overall[,1:3]*(1-output_overall[,1:3])/output_overall[,4])


png("Accuracy_by_method.png", width=5, height=5, units="in", res=300)
coords <- barplot(t(output_overall[,1:3]), beside=TRUE, col=colours, ylim=c(0,1), ylab="Accuracy", las=2)
arrows(x0=coords,y0=t(output_overall[,1:3]), x1=coords, y1=t(output_overall[,1:3]+output_overall_strerr*1.96), length=0)
legend("topright", c("Cyclone", "CycleMix", "Seurat"), fill=colours, bty="n")
dev.off()

labels <- paste0(rownames(time), "\n(",output_overall[,4],")")
png("Time_vs_N.png", width=6, height=6, units="in", res=300)
barplot(t(time), col=colours, ylab="Run Time (s)", beside=TRUE, names=labels)
dev.off()

saveRDS(list(output_overall=output_overall, output_all=output_all, time=time), file="Benchmarking_output.rds")


# By Stage TP, FP, FN
source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")
stage <- "G2M"
TP_tab <- c();
FP_tab <- c();
FN_tab <- c();
ns <- c();
for (d in names(output_all)) {
	n_stage <- output_all[[d]][["cyclone"]][["TPs"]][stage] + output_all[[d]][["cyclone"]][["FNs"]][stage]
	n_other_stage <- output_all[[d]][["cyclone"]]$total - n_stage
	TP_tab <- cbind(TP_tab, c(output_all[[d]][["cyclone"]][["TPs"]][stage]/n_stage, 
				output_all[[d]][["cyclemix"]][["TPs"]][stage]/n_stage, 
				output_all[[d]][["seurat"]][["TPs"]][stage]/n_stage))
	FP_tab <- cbind(FP_tab, c(output_all[[d]][["cyclone"]][["FPs"]][stage]/n_other_stage, 
				output_all[[d]][["cyclemix"]][["FPs"]][stage]/n_other_stage, 
				output_all[[d]][["seurat"]][["FPs"]][stage]/n_other_stage))
	FN_tab <- cbind(FN_tab, c(output_all[[d]][["cyclone"]][["FNs"]][stage]/n_stage, 
				output_all[[d]][["cyclemix"]][["FNs"]][stage]/n_stage,
				output_all[[d]][["seurat"]][["FNs"]][stage]/n_stage))
	ns <- c(ns, n_stage)
}

png(paste("Benchmarking",stage,"bystage.png", sep="_"), width=4, height=8, units="in", res=300)
par(mfrow=c(3,1))

strerr <- t(sqrt(t(TP_tab*(1-TP_tab))/ns))
coords <- barplot(TP_tab, beside=TRUE,names=names(output_all), col=colours, ylim=c(0,1), ylab="TPR", las=2)
arrows(x0=coords,y0=TP_tab, x1=coords, y1=TP_tab+strerr*1.96, length=0)
#legend("topright", c("Cyclone", "CycleMix", "Seurat"), fill=colours, bty="n")

strerr <- t(sqrt(t(FP_tab*(1-FP_tab))/ns))
coords <- barplot(FP_tab, beside=TRUE,names=names(output_all), col=colours, ylab="FPR", las=2, ylim=c(0,1))
arrows(x0=coords,y0=FP_tab, x1=coords, y1=FP_tab+strerr*1.96, length=0)

strerr <- t(sqrt(t(FN_tab*(1-FN_tab))/ns))
coords <- barplot(FN_tab, beside=TRUE,names=names(output_all), col=colours, ylim=c(0,1), ylab="FNR", las=2)
arrows(x0=coords,y0=FN_tab, x1=coords, y1=FN_tab+strerr*1.96, length=0)

dev.off()


