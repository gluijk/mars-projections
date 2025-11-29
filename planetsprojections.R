# Planets projection conversions from equirectangular to orthographic
# www.overfitting.net
# https://www.overfitting.net/2025/11/proyecciones-marcianas-con-r.html

library(tiff)  # save 16-bit TIFF's
library(Rcpp)


# Equirectangular to orthographic projection conversion
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


###########################################################
# MISC PLANETS

planets=c("earth", "moon", "europe", "mars", "venus")

ANCHO=1920
ALTO=ANCHO
for (i in 1:length(planets)) {
    print(paste0("Processing ", planets[i], "..."))
    DEMWHOLE=readTIFF(paste0(planets[i], ".tif"))
    
    for (side in c("L", "R")) {
        if (side=="L") {
            # Left half of DEM
            DEM=DEMWHOLE[,1:(ncol(DEMWHOLE)/2),]
        } else {
            # Right half of DEM
            DEM=DEMWHOLE[,(ncol(DEMWHOLE)/2+1):ncol(DEMWHOLE),]
        }
        
        # Projection conversion
        # imgR, imgG, imgB: full planet texture from -90..+90 longitude and latitude
        out <- equirect_to_orthographic_V2_bilinear(imgR=DEM[,,1], imgG=DEM[,,2], imgB=DEM[,,3],
                                        outW=ANCHO, outH=ALTO, 
                                        Rx=ANCHO/2, Ry=ANCHO/2)
        dim(out) <- c(ALTO, ANCHO, 3)
        writeTIFF(out, paste0("planets_", planets[i], side, ".tif"), bits.per.sample=16)
    }
}




