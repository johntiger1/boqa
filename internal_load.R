boqa.load.data<-function() {
 d<-read.table("/home/johnchen/git_hub/boqa/internal.txt", colClasses=c("integer","integer",rep("numeric",13),"integer"),h=F)
 colnames(d)<-c("run","label","score","marg","marg.ideal", "score.freq","marg.freq", "marg.freq.ideal", "resnik.avg", "resnik.avg.p", "lin.avg", "lin.avg.p", "jc.avg", "jc.avg.p", "mb", "freq");
 return (d);
}
boqa.name<-"/home/johnchen/git_hub/boqa/internal.txt";
boqa.base.name<-"/home/johnchen/git_hub/boqa/internal";
