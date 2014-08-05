#!/usr/bin/R

# Make the mean base quality plot
mnt <- read.table('means.tsv', header=F)
colnames(mnt) <- c('flag', 'meanqual', 'count')
mnt.al.fw <- mnt[mnt$flag == 0,]
mnt.al.rc <- mnt[mnt$flag == 16,]
mnt.unal <- mnt[mnt$flag == 4,]
plot(mnt.al.fw$meanqual, mnt.al.fw$count, typ="o", col="blue", xlab="Mean base quality value", ylab="Count")
points(mnt.al.rc$meanqual, mnt.al.rc$count, typ="o", col="dodgerblue")
points(mnt.unal$meanqual, mnt.unal$count, typ="o", col="red")
legend("topleft", c("Aligned fw strand", "Aligned rc strand", "Unaligned"), col=c("blue", "dodgerblue", "red"), pch=1, lty=1)

# Make the base quality distribution boxplot
qct <- read.table('qual_by_cyc.tsv', header=F)
colnames(qct) <- c('flag', 'cyc', 'qual', 'count')
qct0  <- qct[qct$flag == 0,]
qct4  <- qct[qct$flag == 4,]
qct16 <- qct[qct$flag == 16,]
al0.bx  <- list()
al4.bx  <- list()
al16.bx <- list()
for(i in 0:max(qct$cyc)) {
	ii <- as.character(i)
	al0.bx[[ii]] <- c()
	al4.bx[[ii]] <- c()
	al16.bx[[ii]] <- c()
	for(q in 0:40) {
		qctf0  <- qct0 [qct0$cyc  == i & qct0$qual  == q,]
		qctf4  <- qct4 [qct4$cyc  == i & qct4$qual  == q,]
		qctf16 <- qct16[qct16$cyc == i & qct16$qual == q,]
		if(nrow(qctf0) == 1) {
			al0.bx[[ii]]  <- append(al0.bx[[ii]],  rep(q, round(qctf0$count/10000.0)))
		}
		if(nrow(qctf4) == 1) {
			al4.bx[[ii]]  <- append(al4.bx[[ii]],  rep(q, round(qctf4$count/10000.0)))
		}
		if(nrow(qctf16) == 1) {
			print(round(qctf16$count/10000.0))
			al16.bx[[ii]] <- append(al16.bx[[ii]], rep(q, round(qctf16$count/10000.0)))
		}
	}
}
boxplot(al0.bx)
# rep(qct[11,]$qual, round(qct[11,]$count/10000.0))