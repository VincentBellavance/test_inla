
#---------- SETUP ----------#

library(raster)
library(INLA)
library(ggplot2)

# Variables
species <- "cathartes_aura"
year_start <- 1990 # Première année de données
year_end <- 2020 # Dernière année de données
window_width <- 5 # Largeur de la fenêtre mobile

# Objet spatial propre à chaque espèce
study_extent <- readRDS(paste0("data/",species,"/study_extent.rds"))
obs_all <- readRDS(paste0("data/",species,"/obs.rds"))
rast <- raster::raster(paste0("data/",species,"/rast.gri")) # Raster pour aggréger les données par cellule
qc <- readRDS("data/qc_spacePoly.rds")



#---------- MESH + SPDE ----------#

# Code de François Rousseu
pedge <- 0.03
edge <- min(c(diff(raster::bbox(study_extent)[1,])*pedge,
              diff(raster::bbox(study_extent)[2,])*pedge))
mesh <- INLA::inla.mesh.2d(boundary = study_extent,
                           max.edge = c(edge, edge*5), 
                           min.angle = 21,
                           offset = c(edge, edge*2), 
                           cutoff = edge/2, 
                           crs = raster::crs(study_extent))

spde <- INLA::inla.spde2.pcmatern(mesh=mesh,
                                  alpha=2,
                                  prior.range=c(1000, 0.95),
                                  prior.sigma=c(15, 0.05))



#---------- MODÈLE ET MAP POUR CHAQUE ANNÉE----------#

# Première fenêtre
years <- year_start:(year_start+window_width-1)

# Rouler modèle pour chaque année
for(j in 0:(year_end-years[length(years)])) {

  ## Année du modèle
  year <- mean((years+j))


  #---------- CONVERTIR OBSERVATION EN ESSAI ET SUCCÈS ----------#
  
  ## Filtrer observations pour la fenêtre mobile
  obs <- obs_all[obs_all$year_obs %in% (years+j),]
  
  ## Raster pour stocker succès ou essais par cellule
  succ = rast
  succ[] = 0
  trials = rast
  trials[] = 0
  
  ## Compte de présences et observations totales par cellules
  counts_succ = table(raster::cellFromXY(rast,obs[obs$occurrence == 1,]))
  counts_trials = table(raster::cellFromXY(rast,obs))

  ## Remplir le raster avec les comptes
  succ[as.numeric(names(counts_succ))] = counts_succ
  trials[as.numeric(names(counts_trials))] = counts_trials

  ## Faire un stack avec les deux rasters et convertir en SPDF
  obs <- raster::stack(succ, trials)
  names(obs) <- c("presences", "observations")
  obs <- as(obs, "SpatialPointsDataFrame")

  ## Enlever cellules où il n'y aucun trials
  obs <- obs[obs$observations > 0,]


  #---------- Construire le stack INLA ----------#  

  ## Index matrix
  field <- inla.spde.make.index("field", n.spde=spde$n.spde)

  ## Matrice des coordonnées des observations
	if("x" %in% colnames(as.data.frame(obs_all))){
	  xy <- as.matrix(cbind(obs$x, obs$y))
        } else {
	  xy <- as.matrix(cbind(obs$lon, obs$lat))
	}
	colnames(xy) <- c("x", "y")

  ## A Matrix pour l'estimation
  AEst <- inla.spde.make.A(mesh, loc=xy)
  ## A Matrix pour la prédiction
  APred <- inla.spde.make.A(mesh)

  ## A Matrix dans une liste
  AEstlist <- list(AEst)
  APredlist <- list(APred)

  ## Effets dans une liste (ici seulement variables field)
  effectEst <- list(field)
  effectPred <- list(field)

  ## Constuire le Stack
  ## Pour l'estimation
  StackEst <- inla.stack(data=list(presences = obs$presences,
                                   observations = obs$observations),
                         A = AEstlist,
                         effects = effectEst,
                         tag="est")
  ## Pour la prédiction
  StackPred <- inla.stack(data=list(presences = NA,
                                    observations = NA),
                          A=APredlist,
                          effects=effectPred,
                          tag="pred")
  ## Organiser dans un seul Stack
  Stack <- inla.stack(StackEst,StackPred)


  #---------- Rouler les modèles ----------#

  ## Après la première année, utiliser control.mode pour utiliser le modèle de l'année précédente comme prior
  if(j != 0) {
    model <- inla(presences ~ 0 + f(field, model = spde),
                  data = inla.stack.data(Stack),
                  family="binomial",
                  Ntrials=observations,
                  control.family =list(link="logit"),
                  control.compute=list(waic=TRUE,
                                       openmp.strategy = "huge",
		                                   config = TRUE),
                  control.predictor=list(A=inla.stack.A(Stack),
                                         compute=TRUE, 
                                         link = 1),
                  control.mode = list(theta = previous_model$mode$theta,
                                      restart = TRUE),
                  control.inla=list(int.strategy = "ccd"),
                  verbose = TRUE,
                  debug = TRUE)
  } else {
    model <- inla(presences ~ 0 + f(field, model = spde),
                  data = inla.stack.data(Stack),
                  family="binomial",
                  Ntrials=observations,
                  control.family =list(link="logit"),
                  control.compute=list(waic=TRUE,
                                       openmp.strategy = "huge",
		              config = TRUE),
                  control.predictor=list(A=inla.stack.A(Stack),
                                         compute=TRUE,
                                         link = 1),
                  control.inla=list(int.strategy = "ccd"),
                  verbose = TRUE,
                  debug = TRUE)
  }

  ## Sauvegarder les modèles et les Stacks
  saveRDS(model, paste0("output/",species,"/mod/",year,".rds"))
  saveRDS(Stack, paste0("output/",species,"/stack/",year,".rds"))
  
  
  #---------- Faire les cartes (moyenne et 0.5quant) ----------#
  
  mapBasis <- inla.mesh.projector(mesh,
                                  dims = dim(rast)[2:1],
                                  xlim = c(xmin(study_extent), 
                                           xmax(study_extent)),
                                  ylim = c(ymin(study_extent), 
                                           ymax(study_extent)),
                                  crs = mesh$crs)

  ## Trouver les arêtes du mesh pour y faire les prédiction
  ID <- inla.stack.index(Stack, tag="pred")$data
  
  ## Prédictions avec 0.5quant
  mapPred <- inla.mesh.project(mapBasis, 
                               model$summary.fitted.values[["0.5quant"]][ID])

  ## Convertir en raster
  mapRaster <- raster::raster(t(mapPred[,ncol(mapPred):1]),
                                xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                                ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                                crs = mesh$crs)

  ## Crop et mask pour garder le Qc seulement
  mapRaster <- terra::crop(terra::mask(terra::rast(mapRaster), 
                                       terra::vect(qc)), 
                           terra::vect(qc))

  ## Reconvertir en raster
  map <- raster::raster(mapRaster)
  names(map) <- paste0(year)

  ## Stack des cartes de chaque année ensemble
  if(exists("sdms")) {
    sdms <- raster::stack(sdms, map)
  } else {
    sdms <- raster::stack(map)
  }

  ## Si c'est la dernière année, sauvegarder le stack
  if(j == year_end-years[length(years)]) {
    raster::writeRaster(sdms, 
                        paste0("output/",species,"/maps_pocc"))
  }

  ## Garder le modèle précédent et nettoyer
  previous_model <- model
  rm(model,
     Stack,
     map,
     mapBasis,
     ID,
     mapPred,
     mapRaster)
}



#---------- Calculer le BDI ----------#

# Dataframe avec années, somme des p(occ) et valeur du bdi
bdi <- data.frame(years = 1992:2018,
                  sum_pocc = raster::cellStats(sdms, stat = "sum"),
                  index = c(1, rep(NA, 26)))

# Pas de besoin de faire de log ratio puisqu'on calcul le BDI par espèce
for(i in 2:nrow(bdi)) {
  bdi[i,"index"] <- bdi[i-1,"index"]*(bdi[i,"sum_pocc"]/
                                      bdi[i-1,"sum_pocc"])
}

write.csv(bdi, paste0("output/",species,"/bdi.csv"))



#---------- Faire la figure ----------#

ybreaks <- c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 30) |>
   (\(.) which(. <= max(bdi$index) &
               . >= min(bdi$index)))() |>
     (\(.) seq(min(.)-1, max(.)+1, by = 1))() |>
       (\(.) c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 30)[.])()

p <- ggplot(bdi) +
     geom_line(aes(x = years, 
                   y = index), 
                   lwd = 0.7) +
     geom_hline(yintercept = 1,
                lty = 1,
                col = "grey20",
                lwd = 0.2) +
     labs(y = "Biodiversity Distribution Index", 
          x = "Years") +
     scale_y_continuous(breaks = c(1,ybreaks),
                        trans = "log", 
                        limits = c(min(ybreaks), max(ybreaks))) +
     scale_x_continuous(breaks = c(1992, 2000, 2010, 2018)) +
     theme_bw(base_size = 23)

ggplot2::ggsave(filename = paste0("bdi_",species,".png"),
                plot = p,
                device = "png",
                path = paste0("output/",species),
                width = 15,
                height = 7,
                units = "in",
                bg = "transparent")
