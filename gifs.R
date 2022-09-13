library(raster)

sdms <- stack("progne_subis/maps_pocc.gri")

png(file="progne_subis/img/prosub%02d.png", width = 500, height = 500)
for(i in 1:27) {
  plot(sdms[[i]], main = paste0(c(1992:2018)[i]), cex.main = 2, zlim = c(0,round(max(sdms[,,],na.rm = T), digits = 1)))
}
dev.off()

system("convert -delay 90 progne_subis/img/*.png progne_subis/hirondelle.gif")


sdms <- stack("cathartes_aura/maps_pocc.gri")

png(file="cathartes_aura/img/cataur%02d.png", width = 500, height = 500)
for(i in 1:27) {
  plot(sdms[[i]], main = paste0(c(1992:2018)[i]), cex.main = 2, zlim = c(0,round(max(sdms[,,],na.rm = T), digits = 1)))
}
dev.off()

system("convert -delay 90 cathartes_aura/img/*.png cathartes_aura/urubu.gif")


sdms <- stack("dolichonyx_oryzivorus/maps_pocc.gri")

png(file="dolichonyx_oryzivorus/img/dolory%02d.png", width = 500, height = 500)
for(i in 1:27) {
  plot(sdms[[i]], main = paste0(c(1992:2018)[i]), cex.main = 2, zlim = c(0,round(max(sdms[,,],na.rm = T), digits = 1)))
}
dev.off()

system("convert -delay 90 dolichonyx_oryzivorus/img/*.png dolichonyx_oryzivorus/goglu.gif")