# Natural Earth projection conversions from equirectangular to orthographic
# www.overfitting.net
# https://www.overfitting.net/2025/11/proyecciones-marcianas-con-r.html

library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's
library(Rcpp)


# Equirectangular to orthographic projection conversion (3 versions)

# 1. Single radius
cppFunction('
Rcpp::NumericVector equirect_to_orthographic(
        Rcpp::NumericMatrix imgR,
        Rcpp::NumericMatrix imgG,
        Rcpp::NumericMatrix imgB,
        int outW,
        int outH,
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


# 2. Two radius
cppFunction('
Rcpp::NumericVector equirect_to_orthographic_V2(
        Rcpp::NumericMatrix imgR,
        Rcpp::NumericMatrix imgG,
        Rcpp::NumericMatrix imgB,
        int outW,
        int outH,
        double Rx,
        double Ry)
{
    int texH = imgR.nrow();
    int texW = imgR.ncol();

    Rcpp::NumericVector out(outH * outW * 3);
    double cx = (outW - 1) * 0.5;
    double cy = (outH - 1) * 0.5;

    for (int j = 0; j < outH; j++) {
        for (int i = 0; i < outW; i++) {
            int base = j + i * outH;

            // Normalized orthographic projection coords
            double X = (i - cx) / Rx;
            double Y = (cy - j) / Ry;

            // Outside the elliptical orthographic limb
            if (X*X + Y*Y > 1.0) {
                out[base] = 0.0;
                out[base + outH*outW] = 0.0;
                out[base + 2*outH*outW] = 0.0;
                continue;
            }

            double Z = sqrt(std::max(0.0, 1.0 - X*X - Y*Y));

            double lat = asin(Y);
            double lon = atan2(X, Z);

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


# 3. Two radius + bilinear interpolation
cppFunction('
Rcpp::NumericVector equirect_to_orthographic_V2_bilinear(
        Rcpp::NumericMatrix imgR,
        Rcpp::NumericMatrix imgG,
        Rcpp::NumericMatrix imgB,
        int outW,
        int outH,
        double Rx,
        double Ry)
{
    int texH = imgR.nrow();
    int texW = imgR.ncol();

    Rcpp::NumericVector out(outH * outW * 3);
    double cx = (outW - 1) * 0.5;
    double cy = (outH - 1) * 0.5;

    for (int j = 0; j < outH; j++) {
        for (int i = 0; i < outW; i++) {
            int base = j + i * outH;

            double X = (i - cx) / Rx;
            double Y = (cy - j) / Ry;

            if (X*X + Y*Y > 1.0) {
                out[base] = 0.0;
                out[base + outH*outW] = 0.0;
                out[base + 2*outH*outW] = 0.0;
                continue;
            }

            double Z = sqrt(std::max(0.0, 1.0 - X*X - Y*Y));

            double lat = asin(Y);
            double lon = atan2(X, Z);

            double u = (lon + M_PI/2) / M_PI;
            double v = (M_PI/2 - lat) / M_PI;

            double px = u * (texW - 1);
            double py = v * (texH - 1);

            // Bilinear interpolation indices
            int x0 = std::floor(px);
            int y0 = std::floor(py);
            int x1 = std::min(x0 + 1, texW - 1);
            int y1 = std::min(y0 + 1, texH - 1);

            double dx = px - x0;
            double dy = py - y0;

            // Interpolation weights
            double w00 = (1 - dx) * (1 - dy);
            double w10 = dx * (1 - dy);
            double w01 = (1 - dx) * dy;
            double w11 = dx * dy;

            // R channel
            double R00 = imgR(y0, x0);
            double R10 = imgR(y0, x1);
            double R01 = imgR(y1, x0);
            double R11 = imgR(y1, x1);
            out[base] =
                R00 * w00 + R10 * w10 + R01 * w01 + R11 * w11;

            // G channel
            double G00 = imgG(y0, x0);
            double G10 = imgG(y0, x1);
            double G01 = imgG(y1, x0);
            double G11 = imgG(y1, x1);
            out[base + outH*outW] =
                G00 * w00 + G10 * w10 + G01 * w01 + G11 * w11;

            // B channel
            double B00 = imgB(y0, x0);
            double B10 = imgB(y0, x1);
            double B01 = imgB(y1, x0);
            double B11 = imgB(y1, x1);
            out[base + 2*outH*outW] =
                B00 * w00 + B10 * w10 + B01 * w01 + B11 * w11;
        }
    }
    return out;
}
')


draw_dotted_grid <- function(img, N, dot_step, thickness) {
    DIMX <- dim(img)[1]
    DIMY <- dim(img)[2]
    
    # Positions of lines (excluding borders)
    xs <- round(seq(1, DIMX, length.out = N + 2))[-c(1, N + 2)]
    ys <- round(seq(1, DIMY, length.out = N + 2))[-c(1, N + 2)]
    
    half_t <- floor(thickness / 2)
    
    # Dotted vertical lines
    for (x in xs) {
        for (y in seq(1, DIMY, by = dot_step)) {
            
            x0 <- max(1, x - half_t)
            x1 <- min(DIMX, x + half_t)
            y0 <- max(1, y - half_t)
            y1 <- min(DIMY, y + half_t)
            
            img[x0:x1, y0:y1, ] <- 1
        }
    }
    
    # Dotted horizontal lines
    for (y in ys) {
        for (x in seq(1, DIMX, by = dot_step)) {
            
            x0 <- max(1, x - half_t)
            x1 <- min(DIMX, x + half_t)
            y0 <- max(1, y - half_t)
            y1 <- min(DIMY, y + half_t)
            
            img[x0:x1, y0:y1, ] <- 1
        }
    }
    
    img
}


###########################################################
# VINTAGE EARTH GLOBE


# READ raster
DEMWHOLE=readTIFF("naturalearth.tif")  # read TIFF (not GeoTIFF) file

# Central area of DEM
MIDDLE=ncol(DEMWHOLE)/2
DEM=DEMWHOLE[, (MIDDLE-MIDDLE/2+1):(MIDDLE+MIDDLE/2), ]


# Add grid and shadow
grid <- draw_dotted_grid(DEM*0, 
                         N = 9,          # lines
                         dot_step = 50,  # dot spacing
                         thickness = 10  # pixels thickness
)
writeTIFF(grid, "grid.tif", bits.per.sample=8)
DEM[grid==1]=1

shadow=readTIFF("shadow.tif")
shadow=shadow+0.25
shadow[shadow>1]=1
DEM=DEM*shadow

writeTIFF(DEM, "DEMearth.tif", bits.per.sample=16)


# Projection conversion

# imgR, imgG, imgB: full planet texture from -90..+90 longitude and latitude
ANCHO=1920*4
ALTO=ANCHO
out <- equirect_to_orthographic_V2_bilinear(imgR=DEM[,,1], imgG=DEM[,,2], imgB=DEM[,,3],
                                   outW=ANCHO, outH=ALTO, 
                                   Rx=ANCHO/2, Ry=ALTO/2)
dim(out) <- c(ALTO, ANCHO, 3)
writeTIFF(out, "mapping_earth4.tif", bits.per.sample=16)
