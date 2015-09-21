get.contacts <- function(antlist){
    n=length(antlist)
    contact.list=integer()
    for(i in 1:(n-1)){
        cat("\n", i,"\n")
        for(j in (i+1):n){
            cat(j," ")
            for(k in 1:length(antlist[[i]]$cells)){
                possible=which(antlist[[j]]$cells==antlist[[i]]$cells[k])
                if(length(possible)>0){
                    start.time.i=antlist[[i]]$time[k]
                    end.time.i=antlist[[i]]$time[k]+antlist[[i]]$t.in.cell[k]
                    possible.start=antlist[[j]]$time[possible]
                    possible.end=possible.start+antlist[[j]]$t.in.cell[possible]
                    idx.contact=which(possible.start<end.time.i & possible.end>start.time.i)
                    if(length(idx.contact)>0){
                        contact.list=rbind(contact.list,cbind(i,j,
                            apply(cbind(start.time.i[idx.contact],possible.start[idx.contact]),1,max),
                            apply(cbind(end.time.i[idx.contact],possible.end[idx.contact]),1,min))
                        )
                    }
                }
            }
        }
    }
##    browser()
    contact.list=cbind(contact.list,contact.list[,4]-contact.list[,3])
    colnames(contact.list) <- c("i","j","start","end","length")
    contact.list
}
                        
