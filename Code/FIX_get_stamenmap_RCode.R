# get_stamenmap was not able to download terrain-background tiles.  It was calling them jpeg files when they were actually png files.
# to fix this I modified the get_stamenmap function so that this map type would be considered a png instead of a jpeg and I changed 
# the lapplyF function to lapply.  You need to copy and paste the modified function and it's dependencies into the R Console.  Check here: 
# https://github.com/dkahle/ggmap/blob/master/R/get_stamenmap.R if there are any other pieces missing.

# Abby Lute 8/26/2016

get_stamenmap <- function (bbox = c(left = -95.80204, bottom = 29.38048, right = -94.92313, 
                                    top = 30.14344), zoom = 10, maptype = c("terrain", "terrain-background", 
                                                                            "terrain-labels", "terrain-lines", "toner", "toner-2010", 
                                                                            "toner-2011", "toner-background", "toner-hybrid", "toner-labels", 
                                                                            "toner-lines", "toner-lite", "watercolor"), crop = TRUE, 
                           messaging = FALSE, urlonly = FALSE, color = c("color", "bw"), 
                           force = FALSE, where = tempdir(), ...) 
{
  args <- as.list(match.call(expand.dots = TRUE)[-1])
  argsgiven <- names(args)
  if ("bbox" %in% argsgiven) {
    if (!(is.numeric(bbox) && length(bbox) == 4)) {
      stop("bounding box improperly specified.  see ?get_openstreetmap", 
           call. = F)
    }
  }
  if ("zoom" %in% argsgiven) {
    if (!(is.numeric(zoom) && length(zoom) == 1 && zoom == 
            round(zoom) && zoom >= 0 && zoom <= 18)) {
      stop("scale must be a postive integer 0-18, see ?get_stamenmap.", 
           call. = F)
    }
  }
  if ("messaging" %in% argsgiven) 
    stopifnot(is.logical(messaging))
  if ("urlonly" %in% argsgiven) 
    stopifnot(is.logical(urlonly))
  if ("checkargs" %in% argsgiven) {
    .Deprecated(msg = "checkargs argument deprecated, args are always checked after v2.1.")
  }
  maptype <- match.arg(maptype)
  color <- match.arg(color)
  if (is.null(names(bbox))) 
    names(bbox) <- c("left", "bottom", "right", "top")
  if (maptype %in% c("watercolor")) {
    filetype <- "jpg"
  }
  else {
    filetype <- "png"
  }
  fourCorners <- expand.grid(lon = c(bbox["left"], bbox["right"]), 
                             lat = c(bbox["bottom"], bbox["top"]))
  fourCorners$zoom <- zoom
  row.names(fourCorners) <- c("lowerleft", "lowerright", "upperleft", 
                              "upperright")
  fourCornersTiles <- apply(fourCorners, 1, function(v) LonLat2XY(v[1], 
                                                                  v[2], v[3]))
  xsNeeded <- Reduce(":", sort(unique(as.numeric(sapply(fourCornersTiles, 
                                                        function(df) df$X)))))
  numXTiles <- length(xsNeeded)
  ysNeeded <- Reduce(":", sort(unique(as.numeric(sapply(fourCornersTiles, 
                                                        function(df) df$Y)))))
  numYTiles <- length(ysNeeded)
  tilesNeeded <- expand.grid(x = xsNeeded, y = ysNeeded)
  if (nrow(tilesNeeded) > 40) {
    message(paste0(nrow(tilesNeeded), " tiles needed, this may take a while ", 
                   "(try a smaller zoom)."))
  }
  base_url <- "http://tile.stamen.com/"
  base_url <- paste(base_url, maptype, "/", zoom, sep = "")
  urls <- paste(base_url, apply(tilesNeeded, 1, paste, collapse = "/"), 
                sep = "/")
  urls <- paste(urls, filetype, sep = ".")
  if (messaging) 
    message(length(urls), " tiles required.")
  if (urlonly) 
    return(urls)
  count <- 0
  nTiles <- nrow(tilesNeeded)
  listOfTiles <- lapply(split(tilesNeeded, 1:nrow(tilesNeeded)), 
                        function(v) {
                          v <- as.numeric(v)
                          get_stamenmap_tile(maptype, zoom, v[1], v[2], force = force, 
                                             messaging = messaging)
                        }) # changed lapplyF to lapply 08/26/2016
  map <- stitch(listOfTiles)
  if (!crop) 
    return(map)
  if (crop) {
    mbbox <- attr(map, "bb")
    size <- 256 * c(length(xsNeeded), length(ysNeeded))
    slon <- seq(mbbox$ll.lon, mbbox$ur.lon, length.out = size[1])
    slat <- vector("double", length = 256 * length(ysNeeded))
    for (k in seq_along(ysNeeded)) {
      slat[(k - 1) * 256 + 1:256] <- sapply(as.list(0:255), 
                                            function(y) {
                                              XY2LonLat(X = xsNeeded[1], Y = ysNeeded[k], 
                                                        zoom, x = 0, y = y)$lat
                                            })
    }
    slat <- rev(slat)
    keep_x_ndcs <- which(bbox["left"] <= slon & slon <= bbox["right"])
    keep_y_ndcs <- sort(size[2] - which(bbox["bottom"] <= 
                                          slat & slat <= bbox["top"]))
    croppedmap <- map[keep_y_ndcs, keep_x_ndcs]
  }
  croppedmap <- as.raster(croppedmap)
  class(croppedmap) <- c("ggmap", "raster")
  attr(croppedmap, "bb") <- data.frame(ll.lat = bbox["bottom"], 
                                       ll.lon = bbox["left"], ur.lat = bbox["top"], ur.lon = bbox["right"])
  croppedmap
}



get_stamenmap_tile <- function(maptype, zoom, x, y, force = FALSE, messaging = TRUE, where = tempdir()){
  
  # check arguments
  is.wholenumber <-
    function (x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  
  stopifnot(is.wholenumber(zoom) || !(zoom %in% 1:20))
  stopifnot(is.wholenumber(x) || !(0 <= x && x < 2^zoom))
  stopifnot(is.wholenumber(y) || !(0 <= y && y < 2^zoom))
  
  # format url http://tile.stamen.com/[maptype]/[zoom]/[x]/[y].jpg
  if(maptype %in% c("watercolor")){
    filetype <- "jpg"
  } else {
    filetype <- "png"
  }
  url <- paste0(paste0(c("http://tile.stamen.com", maptype, zoom, x, y), collapse = "/"), ".", filetype)
  
  # lookup in archive
  tile <- file_drawer_get(url)
  if (!is.null(tile) && !force) return(tile)
  
  # grab if not in archive
  tmp <- tempfile()
  download.file(url, destfile = tmp, quiet = !messaging, mode = "wb")
  if(TRUE) message(paste0("Map from URL : ", url))
  
  # read in
  if(maptype %in% c("watercolor")){
    tile <- readJPEG(tmp)
  } else {
    tile <- readPNG(tmp)
  }
  
  
  # convert to colors
  # toner-lines treated differently for alpha
  if(maptype %in% c("toner-hybrid", "toner-labels", "toner-lines",
                    "terrain-labels", "terrain-lines")){
    tile <- t(apply(tile, 1:2, function(x) rgb(x[1], x[2], x[3], x[4])))
  } else {
    tile <- t(apply(tile, 2, rgb))
  }
  
  
  # determine bbox of map. note : not the same as the argument bounding box -
  # the map is only a covering of the bounding box extent the idea is to get
  # the lower left tile and the upper right tile and compute their bounding boxes
  # tiles are referenced by top left of tile, starting at 0,0
  # see http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
  lonlat_upperleft <- XY2LonLat(x, y, zoom)
  lonlat_lowerright <- XY2LonLat(x, y, zoom, 255, 255)
  bbox <- c(
    left = lonlat_upperleft$lon,
    bottom = lonlat_lowerright$lat,
    right = lonlat_lowerright$lon,
    top = lonlat_upperleft$lat
  )
  bb <- data.frame(
    ll.lat = unname(bbox["bottom"]),
    ll.lon = unname(bbox["left"]),
    ur.lat = unname(bbox["top"]),
    ur.lon = unname(bbox["right"])
  )
  
  # format
  class(tile) <- c("ggmap", "raster")
  attr(tile, "bb") <- bb
  
  # store
  file_drawer_set(url, tile)
  
  # return
  tile
}


LonLat2XY <- function(lon_deg, lat_deg, zoom, xpix=256, ypix=256){
  n <- 2^zoom
  X <- ((lon_deg + 180) / 360) * n
  sec <- function(x) 1/cos(x)
  lat_rad <- lat_deg * pi/180
  Y <- (1 - (log(tan(lat_rad) + sec(lat_rad)) / pi)) / 2 * n
  df <- data.frame(
    X = floor(X),
    Y = floor(Y),
    x = xpix*(X - floor(X)),
    y = xpix*(Y - floor(Y))
  )
  row.names(df) <- NULL
  df
}

XY2LonLat <- function(X, Y, zoom, x = 0, y = 0, xpix=255, ypix=255){
  n <- 2^zoom
  lon_deg <- (X+x/xpix) / n * 360.0 - 180.0
  tmp <- tanh( pi * (1 - 2 * (Y+y/ypix) / n))
  ShiftLat <- function(tmp) {
    lat <- 2 * pi * (-1:1) + asin(tmp)
    lat[which(-pi/2 < lat & lat <= pi/2)] * 180/pi
  }
  lat_deg <- ShiftLat(tmp)
  data.frame(lon = lon_deg, lat = lat_deg)
}

library(png)
library(plyr)

#' Manage the ggmap file drawer.
#'
#' To avoid downloading the same maps again and again, ggmap uses a
#' "file drawer" to store the maps you've downloaded. By default,
#' they will be stored in a session specific temporary directory
#' which will be deleted each time you quit R. If you want them to
#' persist across sessions, set the global option
#' "ggmap.file_drawer".
#'
#' @export
#' @keywords internal
#' @examples
#'
#' \dontrun{ if server is unresponsive, R check will fail
#'
#' file_drawer()
#' dir(file_drawer())
#'
#' # The first time you run this, the tiles are downloaded from the server
#' map <- get_stamenmap()
#' # The second time, the local cache is used so it's much faster
#' map <- get_stamenmap()
#'
#' dir(file_drawer())
#'
#' }
#'
file_drawer <- function(...) {
  fd <- getOption("ggmap.file_drawer", file.path(tempdir(), "ggmap"))
  file.path(fd, ...)
}

file_drawer_exists <- function() {
  file.exists(file_drawer())
}

file_drawer_create <- function() {
  if (file_drawer_exists()) return()
  
  dir.create(file_drawer(), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(), file_drawer("index.rds"))
  
  invisible(TRUE)
}

file_drawer_index <- function() {
  file_drawer_create()
  readRDS(file_drawer("index.rds"))
}

file_drawer_set <- function(url, map, name = NULL) {
  if (is.null(name)) {
    name <- paste0(digest::digest(url), '.rds')
  }
  
  index <- file_drawer_index()
  
  if (url %in% names(index)) {
    file.remove(index[[url]])
  }
  index[[url]] <- name
  saveRDS(index, file_drawer("index.rds"))
  saveRDS(map, file_drawer(name))
  
  invisible(TRUE)
}

file_drawer_get <- function(url) {
  index <- file_drawer_index()
  name <- index[[url]]
  
  if (is.null(name)) return(NULL)
  readRDS(file_drawer(name))
}


stitch <- function(tiles){
  
  # trick R CMD check
  ll.lat <- NULL; rm(ll.lat);
  ll.lon <- NULL; rm(ll.lon);
  
  # determine bounding box
  bbs <- ldply(tiles, function(x) attr(x, "bb"))
  
  bigbb <- data.frame(
    ll.lat = min(bbs$ll.lat),
    ll.lon = min(bbs$ll.lon),
    ur.lat = max(bbs$ur.lat),
    ur.lon = max(bbs$ur.lon)
  )
  
  # determine positions of tile in slate (aggregate)
  order <- as.numeric( arrange(bbs, desc(ll.lat), ll.lon)$.id )
  tiles <- tiles[order]
  tiles <- lapply(tiles, as.matrix) # essential for cbind/rbind to work properly!
  
  # split tiles, then squeeze together from top and bottom
  # and then squeeze together from left and right
  nrows <- length( unique(bbs$ll.lat) )
  ncols <- length( unique(bbs$ll.lon) )
  tiles <- split(tiles, rep(1:nrows, each = ncols))
  tiles <- lapply(tiles, function(x) Reduce(cbind, x))
  tiles <- Reduce(rbind, tiles)
  
  tiles <- as.raster(tiles)
  class(tiles) <- c("ggmap", "raster")
  attr(tiles, "bb") <- bigbb
  
  tiles
}
