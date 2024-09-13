# Porthmeus
# 25.ß1.24


require(data.table)

# data dir
dat.dir <- file.path("..","data","Host")

# permutations
n <- 99999
# annotations
subsystems <- fread(file.path(dat.dir,"subsystems.csv"))
subsystems <- subsystems[,unique(subsystem)]

# load enriched subsystems
enriched.subs <- fread(file = file.path(dat.dir, "EnrichedSubsystems.csv"))

st <- "HBMayo"
s.blood <- enriched.subs[set == st & tissue == "blood", unique(subsystem)]
s.biopsy <- enriched.subs[set == st & tissue == "biopsy", unique(subsystem)]
n.blood <- length(s.blood)
n.biopsy <- length(s.biopsy)
n.overlap <- length(intersect(s.blood, s.biopsy))
larger <- rep(NA, n)

for(i in 1:n){
    sm.s.blood <- sample(subsystems, n.blood)
    sm.s.biopsy <- sample(subsystems, n.biopsy)
    sm.n.overlap <- length(intersect(sm.s.blood, sm.s.biopsy))
    larger[i] <- sm.n.overlap >= n.overlap
}
table(larger)
p.val <- sum(larger)/n
