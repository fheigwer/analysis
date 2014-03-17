numbers=c(	"003","004","005","006","007"
			,"008","009","010","011"
			,"012","013","014","015"
			,"016","017","018","019"
			,"020","021","022","023"
			,"024")
letter=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")

for(i in numbers){
	for(j in letter){
		#i="005"
		#j="B"
		wellfolder=paste("/Volumes/New Volume/PathwayData/Christian/2011-04-21_000/Well ",j,i,sep="")
		nucimage=readImage(paste(wellfolder,"/Hoechst_CV - n000000.tif",sep=""))[550:2688,]
		nucimage[which(nucimage>quantile(nucimage,0.999))]=mean(nucimage)
		#nucimage=medianFilter(nucimage, 5, cacheSize=512)
		nucimage=(nucimage-min(nucimage))/(max(nucimage)-min(nucimage))
		

		octimage=readImage(paste(wellfolder,"/FITC_CV - n000000.tif",sep=""))[550:2688,]
		octimage[which(octimage>quantile(octimage,0.999))]=mean(octimage)
		#octimage =medianFilter(octimage, 5, cacheSize=512)
		octimage =(octimage-min(octimage))/(max(octimage)-min(octimage))		
		img = rgbImage(green= octimage, blue= nucimage)
		writeImage(img, file=paste("2011-04-21_000",j,i,".tif",sep="_"), type="tiff", quality = 100, 16 )
	}
}