# Mars projection conversions from equirectangular or Mercator to orthographic
# www.overfitting.net
# https://www.overfitting.net/2025/11/proyecciones-marcianas-con-r.html

library(terra)
library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's
library(Rcpp)


# Equirectangular to orthographic projection conversion
cppFunction('
Rcpp::NumericVector equirect_to_orthographic(
        Rcpp::NumericMatrix imgR,
        Rcpp::NumericMatrix imgG,
        Rcpp::NumericMatrix imgB,
        int outH,
        int outW,
        double R)
{
    int texH = imgR.nrow();
    int texW = imgR.ncol();

    Rcpp::NumericVector out(outH * outW * 3);
    double cx = (outW - 1) * 0.5;
    double cy = (outH - 1) * 0.5;

    for (int j = 0; j < outH; j++) {
        for (int i = 0; i < outW; i++) {
            int base = j + i * outH;

            double X = (i - cx) / R;
            double Y = (cy - j) / R;

            if (X*X + Y*Y > 1.0) {
                out[base] = 0.0;
                out[base + outH*outW] = 0.0;
                out[base + 2*outH*outW] = 0.0;
                continue;
            }

            double Z = sqrt(std::max(0.0, 1.0 - X*X - Y*Y));

            // Convert to spherical coordinates
            double lat = asin(Y);          // latitude
            double lon = atan2(X, Z);      // longitude relative to observer

            // Map to texture coordinates
            double u = (lon + M_PI/2) / M_PI;
            double v = (M_PI/2 - lat) / M_PI;

            double px = u * (texW - 1);
            double py = v * (texH - 1);

            int ix = std::min(std::max(int(px + 0.5), 0), texW-1);
            int iy = std::min(std::max(int(py + 0.5), 0), texH-1);

            out[base] = imgR(iy, ix);
            out[base + outH*outW] = imgG(iy, ix);
            out[base + 2*outH*outW] = imgB(iy, ix);
        }
    }
    return out;
}
')


# Mercator to equirectangular projection conversion
mercator_to_equirect <- function(M,
                                 latmin = -90, latmax = 90,
                                 nlon = dim(M)[2],
                                 nlat = dim(M)[1]) {
    
    # Degree/radian helpers
    deg2rad <- function(d) d * pi / 180
    rad2deg <- function(r) r * 180 / pi
    
    # Mercator forward and inverse transforms
    merc_y <- function(phi_deg) {
        phi <- deg2rad(phi_deg)
        log(tan(pi/4 + phi/2))
    }
    merc_inv <- function(y) {
        rad2deg(2 * atan(exp(y)) - pi/2)
    }
    
    # ---- 1. INPUT GRID: real Mercator y values  ----
    y_min <- merc_y(latmin)
    y_max <- merc_y(latmax)
    
    # Input rows are equally spaced in y_mercator
    in_y <- seq(y_min, y_max, length.out = nlat)
    
    # Corresponding (true) latitudes of input rows
    in_lat <- merc_inv(in_y)
    
    # ---- 2. OUTPUT GRID: equirectangular latitudes ----
    out_lat <- seq(latmin, latmax, length.out = nlat)
    
    # Convert target latitudes to Mercator y
    out_y <- merc_y(out_lat)
    
    # ---- 3. Prepare output structure ----
    is_rgb <- length(dim(M)) == 3
    if (is_rgb) {
        nbands <- dim(M)[3]
        out <- array(NA, dim = c(nlat, nlon, nbands))
    } else {
        out <- matrix(NA, nlat, nlon)
    }
    
    # ---- 4. Nearest neighbour remapping ----
    for (i in 1:nlat) {
        
        # Find nearest input row in Mercator space
        src_i <- which.min(abs(in_y - out_y[i]))
        
        # Copy row
        if (is_rgb) {
            out[i, , ] <- M[src_i, , ]
        } else {
            out[i, ] <- M[src_i, ]
        }
    }
    
    out
}


###########################################################
# 1. NASA/ESA MARS DEM


# READ raster
mars=rast("Mars_MGS_MOLA_DEM_mosaic_global_463m.tif")  # read GeoTIFF file
mars
# dimensions  : 23040, 46080 pixels
# resolution  : 463.0935 m/pixel
# extent      : -10669675, 10669675, -5334838, 5334838 m
# coord. ref. : Equirectangular Mars
# min value   : -8201 
# max value   : 21241 

plot(mars)
RESOLUTION=res(mars)[1]

DEMWHOLE=as.matrix(mars, wide=TRUE)
DEMWHOLE=DEMWHOLE-min(DEMWHOLE)
hist(DEMWHOLE, breaks=800)
MAXHEIGHT=max(DEMWHOLE)

# Left half of DEM
DEM=DEMWHOLE[,1:(ncol(DEMWHOLE)/2)]
writeTIFF(DEM/MAXHEIGHT, "DEM1.tif", bits.per.sample=16)

# Right half of DEM
DEM=DEMWHOLE[,(ncol(DEMWHOLE)/2+1):ncol(DEMWHOLE)]
writeTIFF(DEM/MAXHEIGHT, "DEM2.tif", bits.per.sample=16)


# Projection conversion

# imgR, imgG, imgB: full planet texture from -90..+90 longitude and latitude
ANCHO=1920
ALTO=ANCHO

out <- equirect_to_orthographic(imgR=DEM, imgG=DEM, imgB=DEM,
                                outH=ALTO, outW=ANCHO,
                                R=ANCHO/2)
dim(out) <- c(ALTO, ANCHO, 3)
writeTIFF(out/MAXHEIGHT, "mapping_dem1.tif", bits.per.sample=16)
writeTIFF(out/MAXHEIGHT, "mapping_dem2.tif", bits.per.sample=16)



###########################################################
# 2. PERCIVAL LOWELL MAP

percival=readPNG("mars_map_percival_lowell.png")  # clean map (no borders nor redundancies)
E <- mercator_to_equirect(percival, latmin = -70, latmax = 70)
writePNG(E, "percivalequitectangular.png")


percivalonglat=readPNG("percivalequitectangular180lat.png")  # added +20º/-20º N and S poles

# Left half of DEM
DEM=percivalonglat[,1:(ncol(percivalonglat)/2),]
writeTIFF(DEM, "percivalonglat1.tif", bits.per.sample=16)

# Right half of DEM
DEM=percivalonglat[,(ncol(percivalonglat)/2+1):ncol(percivalonglat),]
writeTIFF(DEM, "percivalonglat2.tif", bits.per.sample=16)


# imgR, imgG, imgB: full planet texture from -90..+90 longitude and latitude
ANCHO=ncol(DEM)
ALTO=ANCHO

# Map DEM
out <- equirect_to_orthographic(imgR=DEM[,,1], imgG=DEM[,,2], imgB=DEM[,,3],
                                outH=ALTO, outW=ANCHO,
                                R=ANCHO/2)
dim(out) <- c(ALTO, ANCHO, 3)
writeTIFF(out, "mapping_percivalonglat1.tif", bits.per.sample=16)
writeTIFF(out, "mapping_percivalonglat2.tif", bits.per.sample=16)


