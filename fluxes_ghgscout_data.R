#####################################################
# Calculate gas fluxes (CO2,CH4) from Picarro scout
#
# Gustaf Granath (adapted after Baird and Green)
#####################################################

# Set constants####

#Standard pressure
st.press=100.00	#kPa

#Standard temperature
st.temp = 273.15	#K

#Volume of 1 mol CO2 and CH4 at STP (they are treated as ideal gases)	
#vol.mol = 22.711	#L mol-1
vol.mol =0.022711	#m3 mol-1

#Molecular mass CO2
mol.mass.co2 = 44.01	#g mol-1
#Molecular mass CH4
mol.mass.ch4 = 16.043	#g mol-1

#Internal area of collar 
coll.area = 0.0661#m2 # for Harbo 0.07793113

#Chamber volume
# now included in meta data
#cham.vol = 0.03 #m3 #harbo 0.00956000

# Get files ####
#library(googledrive)
#dd = as_id("https://drive.google.com/drive/folders/1SvLDE_uhfgUWv_hjj6cuY6_xubr9zsWX?usp=sharing")
#files = drive_ls(path = dd, type = "csv")
#tt = read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", files$id[1]))

# first CO2 and CH4 files
fold <- "C:/Picarro/ghg/"  
filenames = list.files(fold, pattern="*.csv")
filenames # check so no weird csv files are included, like the meta data file!!

#filenames = filenames[1:4] # if one wants a subset and not all files

# second meta data. Remember to not have this file in the same fodler as the data files
#eg meta.dat = read.csv("harbo_scout_picarro_meta.csv")
#meta.dat = read.csv("scout_picarro_meta.csv", dec=",",sep=";")

library(readxl)
meta.dat = read.csv("C:/Picarro/data/metadata.csv")

# Fluxes###
library(ggplot2)
plots = list() # to store plots of each measurement
res=list() # flux results
for (i in 1:length(filenames)) {
  dat <- read.csv(paste(fold, filenames[i], sep=""))
  
  if(any(colnames(dat) == "CO2_dry" | colnames(dat) == "CH4_dry")) {
    match.file = unlist(strsplit(filenames[i], "[_]"))[2]
    match.file = unlist(strsplit(match.file, "Z"))
    
    # Ensure the extracted code is padded with leading zeros to match the format in metadata
    match.file = sprintf("%03d", as.numeric(match.file))
    
    row.id = which(as.character(meta.dat$Code) == match.file)
    if(length(row.id) == 0) { cat("Meta data not found for", match.file, ". Not calculated \n"); next }
    
    print(filenames[i])
    
    # Add pressure and chamber temp data
    t1.press = meta.dat$press[row.id][1]
    t2.press = meta.dat$press[row.id][1]
    
    t1.temp = 273.15 + meta.dat$temp.start[row.id][1]
    t2.temp = 273.15 + meta.dat$temp.end[row.id][1]
    
    # Ensure t1.temp and t2.temp are single values
    if(length(t1.temp) != 1 || length(t2.temp) != 1) {
      cat("Error: t1.temp or t2.temp are not single values\n")
      next
    }
    
    # Temp and pressure change over time
    dat$temp = seq(t1.temp, t2.temp, length.out = nrow(dat))
    dat$press = seq(t1.press, t2.press, length.out = nrow(dat))
    
    if("CO2_dry" %in% colnames(dat)) {
      vol.co2.m3 = (dat$CO2_dry / 1000000) * meta.dat$cham.vol[row.id][1]
      vol.co2.stp.m3 = vol.co2.m3 * (dat$press / dat$temp) * (st.temp / st.press)
      mol.co2 = vol.co2.stp.m3 / vol.mol
      mass.co2 = mol.co2 * mol.mass.co2 * 1000
      min.co2 = min(dat$CO2_dry)
      max.co2 = max(dat$CO2_dry)
      reg = lm(mass.co2 ~ dat$Time_s)
      r2 = summary(lm(mass.co2 ~ dat$Time_s))$r.squared
      p.val = summary(lm(mass.co2 ~ dat$Time_s))$coefficients[2,2]
      flux = (coef(reg)[2] / coll.area) * 86400 # mg CO2 m-2 day-1
      pp = ggplot(data.frame(dat$Time_s, mass.co2), aes_string(x = "dat.Time_s", y = "mass.co2")) +
        geom_point() +
        geom_smooth(method = "lm") +
        ggtitle(paste(as.character(dat$DATE)[1], as.character(dat$TIME)[1], filenames[i], "R2=", round(r2,2),
                      "Flux=", flux, meta.dat$Light_dark[row.id], sep = ","))
      plot(pp)
      plots[[filenames[i]]] <- pp
      res[[filenames[i]]] <- list(filename = filenames[i], Code = match.file, GHG = "co2", flux = flux, r2 = r2, p.value = p.val, startTemp = t1.temp, endTemp = t2.temp,
                                  minCO2 = min.co2, maxCO2 = max.co2)
      rm(pp)
    }
    
    if("CH4_dry" %in% colnames(dat)) {
      vol.ch4.m3 = (dat$CH4_dry / 1000000) * meta.dat$cham.vol[row.id][1]
      vol.ch4.stp.m3 = vol.ch4.m3 * (dat$press / dat$temp) * (st.temp / st.press)
      mol.ch4 = vol.ch4.stp.m3 / vol.mol
      mass.ch4 = mol.ch4 * mol.mass.ch4 * 1000
      min.ch4 = min(dat$CH4_dry)
      max.ch4 = max(dat$CH4_dry) # Use max for consistency
      reg = lm(mass.ch4 ~ dat$Time_s)
      r2 = summary(lm(mass.ch4 ~ dat$Time_s))$r.squared
      p.val = summary(lm(mass.ch4 ~ dat$Time_s))$coefficients[2,2]
      flux = (coef(reg)[2] / coll.area) * 86400 # mg CH4 m-2 day-1
      pp = ggplot(data.frame(dat$Time_s, mass.ch4), aes_string(x = "dat.Time_s", y = "mass.ch4")) +
        geom_point() +
        geom_smooth(method = "lm") +
        ggtitle(paste(as.character(dat$DATE)[1], as.character(dat$TIME)[1], filenames[i], "R2=", round(r2,2),
                      "Flux=", flux, meta.dat$Light_dark[row.id], sep = ","))
      plot(pp)
      plots[[filenames[i]]] <- pp
      res[[filenames[i]]] <- list(filename = filenames[i], Code = match.file, GHG = "ch4", flux = flux, r2 = r2, p.value = p.val, startTemp = t1.temp, endTemp = t2.temp,
                                  minCO2 = min.ch4, maxCO2 = max.ch4)
      rm(pp)
    }
    
  } else {
    cat("No CO2 or CH4 data in file", filenames[i])
    next
  }
}

results = as.data.frame(do.call("rbind",res))
results = tidyr::unnest(results, colnames(results))
results

# save plots
pdf("ghg_time_plots.pdf", width = 12, height = 6)
for (i in 1:3) {
  print(plots[[i]])
}
dev.off()

# merge flux data with id and environmental data
res.dat = merge(results, meta.dat, by="Code")

