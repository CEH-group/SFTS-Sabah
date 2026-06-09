############################# Spatial modelling ################################ 


# load packages
library(INLA)
library(fmesher)
library(sf)


# load data
results <- read.csv("SFTS_data.csv")

# left-join geographic coordiates of household location
results2 = results %>% dplyr::select(-1) %>%
  mutate(houseID2 = .$houseID %>% substr(1, 6))

sh_sfts = st_read("Geolocation//iss_gps.shp")
coords_house = sh_sfts %>% st_coordinates() %>% as_tibble %>% dplyr::select(-3)
coords_house = coords_house %>% 
  mutate(houseID = sh_sfts$name ) 

duplicated_id = coords_house %>% pull(houseID) %>% duplicated


coords_house = coords_house[!duplicated_id,]

mydata_spatial = results2 %>%
  left_join(coords_house, by = c("houseID2"="houseID"))

mydata_spatial_naomit = mydata_spatial %>% na.omit
sf_mydata_inla = st_as_sf(x = mydata_spatial_naomit, coords = c("X","Y"), crs = llCRS)
coords_all = sf_mydata_inla %>% st_coordinates() %>% as_tibble







# filtering Kudat and baggi districts
tmCRS = crs("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
mal_map = st_read("mys_admbnda_adm2_unhcr_20210211.shp")
llCRS = st_crs(mal_map)
my_area = mal_map %>% filter(ADM2_EN == "Kudat") %>% st_transform(tmCRS)

id_banggi = which(coords_all$X>117)
id_kudat = which(coords_all$X<117)

sf_mydata_inla = sf_mydata_inla %>% mutate(area = "kudat")
sf_mydata_inla$area[id_banggi] = "banggi"

# Spatial Modelling for Kudat 

sf_mydata_inla_kudat = sf_mydata_inla %>% filter(area == "kudat")
kudat_coords =sf_mydata_inla_kudat %>% st_transform(tmCRS) %>% st_coordinates(sf_mydata_inla_kudat)
extent_roi = extent(sf_mydata_inla_kudat %>% st_buffer(dist = 1000)%>% st_transform(tmCRS) )

## study region
tmCRS = CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
my_area = mal_map %>% filter(ADM2_EN == "Kudat") %>% st_transform(tmCRS) 

X_max_kudat = 491201
Y_max_kudat = 779822
X_min_kudat = 451049
Y_min_kudat = 731149

ext_kudat <- extent(X_min_kudat, X_max_kudat,
                    Y_min_kudat, Y_max_kudat)

st_kudat = my_area %>% as_Spatial %>%  crop(ext_kudat) %>% st_as_sf
st_kudat_roi = my_area %>% as_Spatial %>%  crop(extent_roi) %>% st_as_sf

## INLA projection matrix
## mesh
fm_hexagon_lattice = function (bnd, edge_len = NULL, buffer_n = 0.49, align = "origin", 
                               meta = FALSE) 
{
  crs <- fm_crs(bnd)
  fm_crs(bnd) <- NA
  bbox <- fm_bbox(bnd)
  if (is.null(edge_len)) {
    edge_len <- diff(bbox[[1]])/250
  }
  if (is.character(align)) {
    align <- match.arg(align, c("origin", "bbox", "centroid"))
    if (align == "bbox") {
      align <- c((bbox[[1]][1] + bbox[[1]][2])/2, (bbox[[2]][1] + 
                                                     bbox[[2]][2])/2)
    }
    else if (align == "centroid") {
      align <- sf::st_centroid(bnd)
    }
    else {
      align <- c(0, 0)
    }
  }
  if (inherits(align, c("sf", "sfc", "sfg"))) {
    align <- sf::st_coordinates(sf::st_centroid(align))
    align <- align[, intersect(colnames(align), c("X", "Y", 
                                                  "Z")), drop = FALSE]
  }
  stopifnot(is.numeric(align))
  origin <- align
  h <- edge_len * sqrt(3)/2
  grid_start <- c(floor((bbox[[1]][1] - origin[1])/edge_len), 
                  floor((bbox[[2]][1] - origin[2])/(2 * h)) * 2L)
  grid_end <- c(ceiling((bbox[[1]][2] - origin[1])/edge_len), 
                ceiling((bbox[[2]][2] - origin[2])/(2 * h)) * 2L)
  grid_n <- grid_end - grid_start + 1L
  if (any(grid_n == 0L)) {
    return()
  }
  x_1_ <- origin[1] + seq(grid_start[1] * edge_len, grid_end[1] * 
                            edge_len, length.out = grid_n[1])
  x_2_ <- origin[1] + seq((grid_start[1] + 0.5) * edge_len, 
                          (grid_end[1] - 0.5) * edge_len, length.out = grid_n[1] - 
                            1L)
  y_1_ <- origin[2] + seq(grid_start[2] * h, grid_end[2] * 
                            h, length.out = (grid_n[2] + 1L)/2L)
  y_2_ <- origin[2] + seq((grid_start[2] + 1) * h, (grid_end[2] - 
                                                      1) * h, length.out = (grid_n[2] + 1L)/2L - 1L)
  x_1 <- rep(x_1_, times = length(y_1_))
  x_2 <- rep(x_2_, times = length(y_2_))
  y_1 <- rep(y_1_, each = length(x_1_))
  y_2 <- rep(y_2_, each = length(x_2_))
  mesh_df <- data.frame(x = c(x_1, x_2), y = c(y_1, y_2))
  lattice_sf <- sf::st_as_sf(mesh_df, coords = c("x", "y"), 
                             crs = crs)
  lattice_sfc <- sf::st_as_sfc(lattice_sf)
  bnd_inner <- sf::st_buffer(bnd, dist = -buffer_n * h)
  fm_crs(bnd_inner) <- crs
  pts_inside <- lengths(sf::st_intersects(lattice_sfc, bnd_inner)) != 
    0
  pts_lattice_sfc <- lattice_sfc[pts_inside]
  if (meta) {
    return(list(lattice = pts_lattice_sfc, edge_len = edge_len, 
                bnd_inner = bnd_inner, grid_n = grid_n, align = origin))
  }
  pts_lattice_sfc
}


mesh_kudat <- fm_mesh_2d(
  loc = fm_hexagon_lattice(bnd = st_kudat, edge_len = 1000 ),
  max.edge = c(1000, 3000)
)

## SPDE
spde_kudat <- inla.spde2.pcmatern(
  mesh = mesh_kudat,
  prior.range = c(1000, 0.5),
  prior.sigma =c(1, 0.5)
  
)

## Index
spde_index_kudat <- inla.spde.make.index("spatial", n.spde = spde_kudat$n.spde)
## A matrix for estimation
A_kudat <- inla.spde.make.A(mesh = mesh_kudat, loc = kudat_coords)

## stack for estimation
results_kudat = c(sf_mydata_inla_kudat$results_bi %>% as.numeric)


stack_kudat_est <- inla.stack(
  data = list(y = results_kudat),
  A = list(A_kudat, 1),
  effects = list(
    list(spatial = spde_index_kudat$spatial),
    list(intercept = rep(1, length(results_kudat)))),
  tag = "est"
)

r_temp = raster(ext = ext_kudat, resolution = c(100, 100), crs = tmCRS)
st_kudat$r =999
r_kudat = st_kudat %>% as_Spatial %>% rasterize(r_temp, field = "r")
r_kudat_roi = r_kudat %>% crop(st_kudat_roi) 
coords_pred_kudat = coordinates(r_kudat_roi)[which(values(r_kudat_roi) ==999),]

# A matrix for prediction
A_kudat_pred <- inla.spde.make.A(mesh = mesh_kudat, loc = coords_pred_kudat)

## stack for prediction
stack_kudat_pred <- inla.stack(
  data = list(y = rep(NA, nrow(coords_pred_kudat))),
  A = list(A_kudat_pred, 1),
  effects = list(
    list(spatial = spde_index_kudat$spatial),
    list(intercept = rep(1, nrow(coords_pred_kudat)))),
  tag = "pred"
)

# combined stack
stack_kudat = inla.stack(stack_kudat_est, stack_kudat_pred)

## INLA
formula_null = y ~ -1 + intercept + f(spatial, model = spde_kudat)

model_null <- inla(
  formula_null,
  data = inla.stack.data(stack_kudat),
  family = "binomial",
  control.predictor = list(A = inla.stack.A(stack_kudat), compute = TRUE),
  control.compute = list(dic = TRUE, 
                         waic = TRUE,
                         cpo = TRUE,
                         return.marginals.predictor = TRUE, 
                         config = TRUE), 
  control.inla = list(lincomb.derived.correlation.matrix = TRUE),
  control.fixed = list(correlation.matrix = TRUE))


summary(model_null)

## Assign predicted values in raster
pred_id = inla.stack.index(stack_kudat, "pred")$data 
pred_prob = model_null$summary.linear.predictor[pred_id,"mean"] %>% plogis 
pred_kudat = r_kudat_roi
values(pred_kudat)[which(values(pred_kudat) ==999)] = pred_prob


plot(pred_kudat)



# Spatial Modelling for Banggi
sf_mydata_inla_banggi = sf_mydata_inla %>% filter(area == "banggi")
banggi_coords = sf_mydata_inla_banggi %>% st_transform(tmCRS) %>% st_coordinates

extent_roi = sf_mydata_inla_banggi %>% st_buffer(dist = 1000) %>% st_transform(tmCRS) %>% extent

X_max_banggi = 546691
Y_max_banggi = 821149

X_min_banggi = 501639
Y_min_banggi = 778863

ext_banggi <- extent(X_min_banggi, X_max_banggi,
                    Y_min_banggi, Y_max_banggi)

st_banggi = my_area %>% as_Spatial %>%  crop(ext_banggi) %>% st_as_sf



## Mesh
mesh_banggi <- fm_mesh_2d(
  loc = fm_hexagon_lattice(bnd = st_banggi, edge_len = 1000 ),
  max.edge = c(1500, 2000), 
  offset = -0.2
)

## SPDE
spde_banggi <- inla.spde2.pcmatern(
  mesh = mesh_bangi,
  prior.range = c(1000, 0.5),
  prior.sigma =c(1, 0.5)
  
)

## index
spde_index_banggi <- inla.spde.make.index("spatial", n.spde = spde_banggi$n.spde)
## A matrix for estimation
A_banggi <- inla.spde.make.A(mesh = mesh_bangi, loc = banggi_coords)

## stack for estimation

results_banggi = c(sf_mydata_inla_banggi$results_bi %>% as.numeric)


r_temp = raster(ext = ext_banggi, resolution = c(100, 100), crs = tmCRS)
st_banggi$r = 999
r_banggi = st_banggi %>% as_Spatial %>% rasterize(r_temp, field = "r")
r_banggi_roi = r_bangi %>% crop(extent_roi)
coords_pred_banggi = coordinates(r_bangi_roi)[which(values(r_banggi_roi) ==999),]

stack_bangi_est <- inla.stack(
  data = list(y = results_bangi),
  A = list(A_banggi, 1),
  effects = list(
    list(spatial = spde_index_bangi$spatial),
    list(intercept = rep(1, length(results_bangi)))),
  tag = "est"
)

## A matrix for prediction
A_banggi_pred <- inla.spde.make.A(mesh = mesh_bangi, loc = coords_pred_banggi)

## stack for prediction
stack_bangi_pred <- inla.stack(
  data = list(y = NA),
  A = list(A_banggi_pred, 1),
  effects = list(
    list(spatial = spde_index_bangi$spatial),
    list(intercept = rep(1, nrow(coords_pred_bangi)))),
  tag = "pred"
)

## Combined stack
stack_bangi = inla.stack(stack_bangi_est, stack_bangi_pred)

# INLA
formula_null = y ~ -1 + intercept + f(spatial, model = spde_banggi)

model_null_bangi <- inla(
  formula_null,
  data = inla.stack.data(stack_bangi),
  family = "binomial",
  control.predictor = list(A = inla.stack.A(stack_bangi), compute = TRUE),
  control.compute = list(dic = TRUE, 
                         waic = TRUE,
                         cpo = TRUE,
                         return.marginals.predictor = TRUE, 
                         config = TRUE), 
  control.inla = list(lincomb.derived.correlation.matrix = TRUE),
  control.fixed = list(correlation.matrix = TRUE))

summary(model_null_bangi)


## Assign predicted values in raster
pred_id = inla.stack.index(stack_bangi, "pred")$data 
pred_prob = model_null_bangi$summary.linear.predictor[pred_id,"mean"] %>% plogis 
pred_bangi= r_bangi_roi
values(pred_bangi)[which(values(r_bangi_roi) ==999)] = pred_prob
plot(pred_bangi)




