# Basic HDR tone mapping algorithm using contrast (sigmoid) curves
# www.overfitting.net
# https://www.overfitting.net/

library(tiff)  # save 16-bit TIFF's
library(terra)
library(png)  # save 8-bit PNG's
library(Rcpp)


cppFunction('
Rcpp::NumericVector ortho_project_fulltexture(
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


# Hillshade calculation
hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction. It can be defined in three ways:
    #   a) 1D value indicating light source Azimuth in degrees (0-360)
    #      (0=North, 90=East, 180=South, 270=West)
    #   b) 2D vector indicating light source (X,Y) coordinates
    #   c) 3D vector indicating light source (X,Y,Z) coordinates:
    #      (X=South, Y=East, Z=Up)
    #      dlight=c(0, 2, 3)  # sunrise
    #      dlight=c(0, 0, 1)  # midday
    #      dlight=c(0,-2, 3)  # sunset
    #   NOTE: both in a) and b) a 45º Elevation angle is applied
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    # Deal with lighting direction
    if (length(dlight)==1) dlight=c(-cos(dlight*pi/180), sin(dlight*pi/180))
    if (length(dlight)==2) dlight=c(dlight, (dlight[1]^2+dlight[2]^2)^0.5)
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}


# It is NOT faster than the original R version
cppFunction('
Rcpp::NumericMatrix hillshade_cpp(Rcpp::NumericMatrix DEM, double dx=25.0,
                                 Rcpp::NumericVector dlight = Rcpp::NumericVector::create(0,2,3),
                                 double gamma=1.0)
{
    int DIMY = DEM.nrow();
    int DIMX = DEM.ncol();

    // Handle lighting vector
    if (dlight.size() == 1) {
        double az = dlight[0] * M_PI / 180.0;
        dlight = Rcpp::NumericVector::create(-cos(az), sin(az));
    }
    if (dlight.size() == 2) {
        double z = std::sqrt(dlight[0]*dlight[0] + dlight[1]*dlight[1]);
        dlight.push_back(z);
    }
    double dlightM = std::sqrt(dlight[0]*dlight[0] + dlight[1]*dlight[1] + dlight[2]*dlight[2]);

    // Initialize hillshade matrix
    Rcpp::NumericMatrix hillshade(DIMY, DIMX);

    // Temporary matrix for central difference slopes
    Rcpp::NumericMatrix nx(DIMY-2, DIMX-2);
    Rcpp::NumericMatrix ny(DIMY-2, DIMX-2);
    Rcpp::NumericMatrix nz(DIMY-2, DIMX-2);

    // Compute slope vectors
    for(int j=0; j<DIMY-2; j++){
        for(int i=0; i<DIMX-2; i++){
            nx(j,i) = 2*dx*(DEM(j,i+1) - DEM(j+2,i+1));
            ny(j,i) = 2*dx*(DEM(j+1,i) - DEM(j+1,i+2));
            nz(j,i) = 4*dx*dx;
        }
    }

    // Dot product and magnitude
    for(int j=0; j<DIMY-2; j++){
        for(int i=0; i<DIMX-2; i++){
            double nM = std::sqrt(nx(j,i)*nx(j,i) + ny(j,i)*ny(j,i) + nz(j,i)*nz(j,i));
            double dn = dlight[0]*nx(j,i) + dlight[1]*ny(j,i) + dlight[2]*nz(j,i);
            double value = dn / (dlightM * nM);
            if(value < 0) value = 0.0;
            hillshade(j+1, i+1) = std::pow(value, 1.0/gamma);
        }
    }

    // Fill borders (copy from second row/column)
    for(int i=0; i<DIMX; i++){
        hillshade(0,i) = hillshade(1,i);
        hillshade(DIMY-1,i) = hillshade(DIMY-2,i);
    }
    for(int j=0; j<DIMY; j++){
        hillshade(j,0) = hillshade(j,1);
        hillshade(j,DIMX-1) = hillshade(j,DIMX-2);
    }

    return hillshade;
}
')


mercator_to_equirect <- function(M,
                                 latmin = -70, latmax = 70,
                                 nlon = dim(M)[2],
                                 nlat = dim(M)[1]) {
    
    # Degree↔radian helpers
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
    
    # Input rows are equally spaced in y_mercator:
    in_y <- seq(y_min, y_max, length.out = nlat)
    
    # Corresponding (true) latitudes of input rows:
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
    
    # ---- 4. Nearest-neighbour remapping ----
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
# 1. READ DEM 

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
writeTIFF(DEM/MAXHEIGHT, "DEM1bis.tif", bits.per.sample=16)

# Right half of DEM
DEM=DEMWHOLE[,(ncol(DEMWHOLE)/2+1):ncol(DEMWHOLE)]
writeTIFF(DEM/MAXHEIGHT, "DEM2.tif", bits.per.sample=16)


###########################################################
# 2. CALCULATE HILLSHADE

equatormars=21297  # 21.297 kilometers
halfmars=equatormars/2  # km represented in the projection
RESOLUTION=halfmars/ncol(DEM)  # resolution in km: 0.4621745 km

hillshade=hillshademap(DEM/10, dx=RESOLUTION*10, dlight=c(0,1,1))  # (X,Y,Z) coords
writeTIFF(hillshade, "hillshade2.tif", bits.per.sample=16)


###########################################################
# 3. PROJECTION

# imgR, imgG, imgB: full planet texture from -90..+90 longitude and latitude
ANCHO=1920
ALTO=ANCHO

# Map DEM
out <- ortho_project_fulltexture(imgR=DEM, imgG=DEM, imgB=DEM,
                                 outH=ALTO, outW=ANCHO,
                                 R=ANCHO/2)
dim(out) <- c(ALTO, ANCHO, 3)
writeTIFF(out/MAXHEIGHT, "mapping_dem2.tif", bits.per.sample=16)

# Map hillshade
hillshade=readTIFF("hillshade2.tif")
out <- ortho_project_fulltexture(imgR=hillshade, imgG=hillshade, imgB=hillshade,
                                 outH=ALTO, outW=ANCHO,
                                 R=ANCHO/2)
dim(out) <- c(ALTO, ANCHO, 3)
writeTIFF(out, "mapping_hillshade2.tif", bits.per.sample=16)



percival=readPNG("mars_map_percival_lowell.png")
E <- mercator_to_equirect(percival, latmin = -70, latmax = 70)
writePNG(E, "percivalEquitectangular.png")


percivalonglat=readPNG("percivalequitectangular180lat.png")
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
out <- ortho_project_fulltexture(imgR=DEM[,,1], imgG=DEM[,,2], imgB=DEM[,,3],
                                 outH=ALTO, outW=ANCHO,
                                 R=ANCHO/2)
dim(out) <- c(ALTO, ANCHO, 3)
writeTIFF(out, "mapping_percivalonglat1.tif", bits.per.sample=16)


