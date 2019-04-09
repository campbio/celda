# Adapted originally from the very excellent pheatmap package
# (https://cran.r-project.org/web/packages/pheatmap/index.html)


.lo <- function(rown,
    coln,
    nrow,
    ncol,
    cellHeight = NA,
    cellWidth = NA,
    treeHeightCol,
    treeHeightRow,
    legend,
    annotationRow,
    annotationCol,
    annotationColors,
    annotationLegend,
    annotationNamesRow,
    annotationNamesCol,
    main,
    fontSize,
    fontSizeRow,
    fontSizeCol,
    gapsRow,
    gapsCol,
    ...) {
    # Get height of colnames and length of rownames
    if (!is.null(coln[1]) |
            (!.is.na2(annotationRow) & annotationNamesRow)) {
        if (!is.null(coln[1])) {
            t <- coln
        } else {
            t <- ""
        }
        tw <- strwidth(t, units = "in", cex = fontSizeCol / fontSize)
        if (annotationNamesRow) {
            t <- c(t, colnames(annotationRow))
            tw <- c(tw, strwidth(colnames(annotationRow), units = "in"))
        }
        longestColn <- which.max(tw)
        gp <- list(fontSize = ifelse(longestColn <= length(coln),
            fontSizeCol,
            fontSize), ...)
        colnHeight <- unit(1,
            "grobheight",
            textGrob(t[longestColn],
                rot = 90,
                gp = do.call(gpar, gp))) +
            unit(10, "bigpts")
    } else {
        colnHeight <- unit(5, "bigpts")
    }

<<<<<<< HEAD
    if (!is.null(rown[1])) {
        t <- rown
        tw <- strwidth(t, units = "in", cex = fontSizeRow / fontSize)
        if (annotationNamesCol) {
            t <- c(t, colnames(annotationCol))
            tw <- c(tw, strwidth(colnames(annotationCol), units = "in"))
=======
    if(!is.null(rown[1])){
        t = rown
        tw = strwidth(t, units = 'in', cex = fontsize_row / fontsize)
        if(annotation_names_col){
            t = c(t, colnames(annotation_col))
            tw = c(tw, strwidth(colnames(annotation_col), units = 'in'))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
        }
        longestRown <- which.max(tw)
        gp <- list(fontSize = ifelse(longestRown <= length(rown),
            fontSizeRow,
            fontSize), ...)
        rownWidth <- unit(1,
            "grobwidth",
            textGrob(t[longestRown],
                gp = do.call(gpar, gp))) +
            unit(10, "bigpts")
    } else {
        rownWidth <- unit(5, "bigpts")
    }

<<<<<<< HEAD
    gp <- list(fontSize = fontSize, ...)
=======
    gp = list(fontsize = fontsize, ...)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
    # Legend position
    if (!.is.na2(legend)) {
        longestBreak <- which.max(nchar(names(legend)))
        longestBreak <- unit(1.1,
                "grobwidth",
                textGrob(as.character(names(legend))[longestBreak],
                    gp = do.call(gpar, gp)))
        titleLength <- unit(1.1,
            "grobwidth",
            textGrob("Scale",
                gp = gpar(fontface = "bold",
                    ...)))
        legendWidth <- unit(12, "bigpts") + longestBreak * 1.2
        legendWidth <- max(titleLength, legendWidth)
    } else {
        legendWidth <- unit(0, "bigpts")
    }

    # Set main title height
    if (is.na(main)) {
        mainHeight <- unit(0, "npc")
    } else {
        mainHeight <- unit(1.5,
            "grobheight",
            textGrob(main,
                gp = gpar(fontSize = 1.3 * fontSize,
                    ...)))
    }

    # Column annotations
<<<<<<< HEAD
    textheight <- unit(fontSize, "bigpts")

    if (!.is.na2(annotationCol)) {
        # Column annotation height
        annotColHeight <-
            ncol(annotationCol) *
            (textheight + unit(2, "bigpts")) +
            unit(2, "bigpts")

        # Width of the correponding legend
        t <- c(as.vector(as.matrix(annotationCol)), colnames(annotationCol))
        annotColLegendWidth <- unit(1.2,
            "grobwidth",
            textGrob(t[which.max(nchar(t))],
                gp = gpar(...))) +
            unit(12, "bigpts")
        if (!annotationLegend) {
            annotColLegendWidth <- unit(0, "npc")
=======
    textheight = unit(fontsize, "bigpts")

    if(!is.na2(annotation_col)){
        # Column annotation height
        annot_col_height = ncol(annotation_col) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

        # Width of the correponding legend
        t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
        annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        if(!annotation_legend){
            annot_col_legend_width = unit(0, "npc")
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
        }
    } else {
        annotColHeight <- unit(0, "bigpts")
        annotColLegendWidth <- unit(0, "bigpts")
    }
<<<<<<< HEAD

    # Row annotations
    if (!.is.na2(annotationRow)) {
        # Row annotation width
        annotRowWidth <- ncol(annotationRow) *
            (textheight + unit(2, "bigpts")) +
            unit(2, "bigpts")

        # Width of the correponding legend
        t <- c(as.vector(as.matrix(annotationRow)),
            colnames(annotationRow))
        annotRowLegendWidth <- unit(1.2,
            "grobwidth",
            textGrob(t[which.max(nchar(t))],
                gp = gpar(...))) +
            unit(12,
                "bigpts")

        if (!annotationLegend) {
            annotRowLegendWidth <- unit(0, "npc")
=======
    else{
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }

    # Row annotations
    if(!is.na2(annotation_row)){
        # Row annotation width
        annot_row_width = ncol(annotation_row) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

        # Width of the correponding legend
        t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
        annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        if(!annotation_legend){
            annot_row_legend_width = unit(0, "npc")
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
        }
    } else {
        annotRowWidth <- unit(0, "bigpts")
        annotRowLegendWidth <- unit(0, "bigpts")
    }
<<<<<<< HEAD

    annotLegendWidth <- max(annotRowLegendWidth, annotColLegendWidth)

    # Tree height
    treeHeightCol <- unit(treeHeightCol, "bigpts") + unit(5, "bigpts")
    treeHeightRow <- unit(treeHeightRow, "bigpts") + unit(5, "bigpts")

    # Set cell sizes
    if (is.na(cellWidth)) {
        matWidth <- unit(1, "npc") -
            rownWidth -
            legendWidth -
            treeHeightRow -
            annotRowWidth -
            annotLegendWidth
    } else {
        matWidth <- unit(cellWidth * ncol, "bigpts") +
            length(gapsCol) *
            unit(0, "bigpts")
    }

    if (is.na(cellHeight)) {
        matHeight <- unit(1, "npc") -
            mainHeight -
            colnHeight -
            treeHeightCol -
            annotColHeight
    } else {
        matHeight <- unit(cellHeight * nrow, "bigpts") +
            length(gapsRow) *
            unit(0, "bigpts")
    }

    # Produce gtable
    gt <- gtable(widths = unit.c(treeHeightRow,
            annotRowWidth,
            matWidth,
            rownWidth,
            legendWidth,
            annotLegendWidth),
        heights = unit.c(mainHeight,
            treeHeightCol,
            annotColHeight,
            matHeight,
            colnHeight),
        vp = viewport(gp = do.call(gpar, gp)))

    cw <- convertWidth(matWidth -
        (length(gapsCol) * unit(0, "bigpts")),
        "bigpts", valueOnly = TRUE) / ncol
    ch <- convertHeight(matHeight -
        (length(gapsRow) * unit(0, "bigpts")),
        "bigpts", valueOnly = TRUE) / nrow

    # Return minimal cell dimension in bigpts to decide if borders are drawn
    mindim <- min(cw, ch)

    res <- list(gt = gt, mindim = mindim)
=======
    else{
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }

    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)

    # Tree height
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts")

    # Set cell sizes
    if(is.na(cellwidth)){
        mat_width = unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_row_width - annot_legend_width
    }
    else{
        mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(0, "bigpts")
    }

    if(is.na(cellheight)){
        mat_height = unit(1, "npc") - main_height - coln_height - treeheight_col - annot_col_height
    }
    else{
        mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(0, "bigpts")
    }

    # Produce gtable
    gt = gtable(widths = unit.c(treeheight_row, annot_row_width, mat_width, rown_width, legend_width, annot_legend_width), heights = unit.c(main_height, treeheight_col, annot_col_height, mat_height, coln_height), vp = viewport(gp = do.call(gpar, gp)))

    cw = convertWidth(mat_width - (length(gaps_col) * unit(0, "bigpts")), "bigpts", valueOnly = TRUE) / ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(0, "bigpts")), "bigpts", valueOnly = TRUE) / nrow

    # Return minimal cell dimension in bigpts to decide if borders are drawn
    mindim = min(cw, ch)

    res = list(gt = gt, mindim = mindim)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

.findCoordinates <- function(n, gaps, m = seq(1,n)) {
    if (length(gaps) == 0) {
        return(list(
            coord = unit(m / n, "npc"),
            size = unit(1 / n, "npc")))
    }

<<<<<<< HEAD
    if (max(gaps) > n) {
        stop("Gaps do not match with matrix size")
    }

    size <- (1 / n) *
        (unit(1, "npc") - length(gaps) * unit("0", "bigpts"))

    gaps2 <- base::apply(sapply(gaps,
        function(gap, x) {x > gap},
        m), 1, sum)
    coord <- m * size + (gaps2 * unit("0", "bigpts"))
=======
    if(max(gaps) > n){
        stop("Gaps do not match with matrix size")
    }

    size = (1 / n) * (unit(1, "npc") - length(gaps) * unit("0", "bigpts"))

    gaps2 = base::apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum)
    coord = m * size + (gaps2 * unit("0", "bigpts"))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(list(coord = coord, size = size))
}

.drawDendrogram <- function(hc, gaps, horizontal = TRUE) {
    h <- hc$height / max(hc$height) / 1.05
    m <- hc$merge
    o <- hc$order
    n <- length(o)

<<<<<<< HEAD
    m[m > 0] <- n + m[m > 0]
    m[m < 0] <- abs(m[m < 0])

    dist <- matrix(0,
        nrow = 2 * n - 1,
        ncol = 2,
        dimnames = list(NULL, c("x", "y")))
    dist[seq(1,n), 1] <- 1 / n / 2 + (1 / n) *
        (match(seq(1,n), o) - 1)
=======
    m[m > 0] = n + m[m > 0]
    m[m < 0] = abs(m[m < 0])

    dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y")))
    dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    for (i in seq(1,nrow(m))) {
        dist[n + i, 1] <- (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
        dist[n + i, 2] <- h[i]
    }

<<<<<<< HEAD
    drawConnection <- function(x1, x2, y1, y2, y) {
        res <- list(x = c(x1, x1, x2, x2),
            y = c(y1, y, y, y2))
=======
    draw_connection = function(x1, x2, y1, y2, y){
        res = list(
            x = c(x1, x1, x2, x2),
            y = c(y1, y, y, y2)
        )
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

        return(res)
    }

<<<<<<< HEAD
    x <- rep(NA, nrow(m) * 4)
    y <- rep(NA, nrow(m) * 4)
    id <- rep(1:nrow(m), rep(4, nrow(m)))

    for (i in seq(1,nrow(m))) {
        c <- drawConnection(dist[m[i, 1], 1],
            dist[m[i, 2], 1],
            dist[m[i, 1], 2],
            dist[m[i, 2], 2],
            h[i])
        k <- (i - 1) * 4 + 1
        x[seq(k, k + 3)] <- c$x
        y[seq(k, k + 3)] <- c$y
    }

    x <- .findCoordinates(n, gaps, x * n)$coord
    y <- unit(y, "npc")

    if (!horizontal) {
        a <- x
        x <- unit(1, "npc") - y
        y <- unit(1, "npc") - a
    }
    res <- polylineGrob(x = x, y = y, id = id)
=======
    x = rep(NA, nrow(m) * 4)
    y = rep(NA, nrow(m) * 4)
    id = rep(1:nrow(m), rep(4, nrow(m)))

    for(i in 1:nrow(m)){
        c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
        k = (i - 1) * 4 + 1
        x[k : (k + 3)] = c$x
        y[k : (k + 3)] = c$y
    }

    x = find_coordinates(n, gaps, x * n)$coord
    y = unit(y, "npc")

    if(!horizontal){
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
    }
    res = polylineGrob(x = x, y = y, id = id)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

<<<<<<< HEAD
.drawMatrix <- function(matrix,
    borderColor,
    gapsRows,
    gapsCols,
    fmat,
    fontSizeNumber,
    numberColor) {

    n <- nrow(matrix)
    m <- ncol(matrix)

    coordX <- .findCoordinates(m, gapsCols)
    coordY <- .findCoordinates(n, gapsRows)

    x <- coordX$coord -
        0.5 * coordX$size
    y <- unit(1, "npc") -
        (coordY$coord - 0.5 * coordY$size)

    coord <- expand.grid(y = y, x = x)

    res <- gList()

    res[["rect"]] <- rectGrob(x = coord$x,
        y = coord$y,
        width = coordX$size,
        height = coordY$size,
        gp = gpar(fill = matrix, col = borderColor))

    if (attr(fmat, "draw")) {
        res[["text"]] <- textGrob(x = coord$x,
            y = coord$y,
            label = fmat,
            gp = gpar(col = numberColor, fontSize = fontSizeNumber))
    }

    res <- gTree(children = res)
=======
draw_matrix = function(matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color){
    n = nrow(matrix)
    m = ncol(matrix)

    coord_x = find_coordinates(m, gaps_cols)
    coord_y = find_coordinates(n, gaps_rows)

    x = coord_x$coord - 0.5 * coord_x$size
    y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)

    coord = expand.grid(y = y, x = x)

    res = gList()

    res[["rect"]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = coord_y$size, gp = gpar(fill = matrix, col = border_color))

    if(attr(fmat, "draw")){
        res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
    }

    res = gTree(children = res)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

<<<<<<< HEAD
.drawColnames <- function(coln, gaps, ...) {
    coord <- .findCoordinates(length(coln), gaps)
    x <- coord$coord - 0.5 * coord$size

    res <- textGrob(coln,
        x = x,
        y = unit(1, "npc") -
            unit(3, "bigpts"),
        vjust = 0.5,
        hjust = 0,
        rot = 270,
        gp = gpar(...))
=======
draw_colnames = function(coln, gaps, ...){
    coord = find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size

    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

<<<<<<< HEAD
.drawRownames <- function(rown, gaps, ...) {
    coord <- .findCoordinates(length(rown), gaps)
    y <- unit(1, "npc") - (coord$coord - 0.5 * coord$size)

    res <- textGrob(rown,
            x = unit(3, "bigpts"),
            y = y,
            vjust = 0.5,
            hjust = 0,
            gp = gpar(...))
=======
draw_rownames = function(rown, gaps, ...){
    coord = find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)

    res = textGrob(rown, x = unit(3, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

<<<<<<< HEAD
.drawLegend <- function(color, breaks, legend, ...) {
    height <- min(unit(1, "npc"), unit(150, "bigpts"))

    legendPos <- (legend - min(breaks)) / (max(breaks) - min(breaks))
    legendPos <- height * legendPos + (unit(1, "npc") - height)

    breaks <- (breaks - min(breaks)) / (max(breaks) - min(breaks))
    breaks <- height * breaks + (unit(1, "npc") - height)

    h <- breaks[-1] - breaks[-length(breaks)]

    rect <- rectGrob(x = 0,
        y = breaks[-length(breaks)],
        width = unit(10, "bigpts"),
        height = h,
        hjust = 0,
        vjust = 0,
        gp = gpar(fill = color, col = "#FFFFFF00"))

    text <- textGrob(names(legend),
        x = unit(14, "bigpts"),
        y = legendPos,
        hjust = 0,
        gp = gpar(...))

    res <- grobTree(rect, text)
=======
draw_legend = function(color, breaks, legend, ...){
    height = min(unit(1, "npc"), unit(150, "bigpts"))

    legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
    legend_pos = height * legend_pos + (unit(1, "npc") - height)

    breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
    breaks = height * breaks + (unit(1, "npc") - height)

    h = breaks[-1] - breaks[-length(breaks)]

    rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
    text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))

    res = grobTree(rect, text)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

<<<<<<< HEAD
.convertAnnotations <- function(annotation, annotationColors) {
    new <- annotation
    for (i in 1:ncol(annotation)) {
        a <- annotation[, i]
        b <- annotationColors[[colnames(annotation)[i]]]
        if (is.character(a) | is.factor(a)) {
            a <- as.character(a)

            if (length(setdiff(setdiff(a, NA), names(b))) > 0) {
                stop(sprintf("Factor levels on variable %s do not match
                        with annotationColors",
                    colnames(annotation)[i]))
=======
convert_annotations = function(annotation, annotation_colors){
    new = annotation
    for(i in 1:ncol(annotation)){
        a = annotation[, i]
        b = annotation_colors[[colnames(annotation)[i]]]
        if(is.character(a) | is.factor(a)){
            a = as.character(a)

            if(length(setdiff(setdiff(a, NA), names(b))) > 0){
                stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
            }
            new[, i] <- b[a]
        } else {
            a <- cut(a, breaks = 100)
            new[, i] <- colorRampPalette(b)(100)[a]
        }
    }
    return(as.matrix(new))
}

<<<<<<< HEAD
.drawAnnotations <- function(convertedAnnotations,
    borderColor,
    gaps,
    fontSize,
    horizontal) {

    n <- ncol(convertedAnnotations)
    m <- nrow(convertedAnnotations)

    coordX <- .findCoordinates(m, gaps)

    x <- coordX$coord - 0.5 * coordX$size

    # y = cumsum(rep(fontSize, n)) - 4 + cumsum(rep(2, n))
    y <- cumsum(rep(fontSize, n)) +
        cumsum(rep(2, n)) -
        fontSize / 2 + 1
    y <- unit(y, "bigpts")

    if (horizontal) {
        coord <- expand.grid(x = x, y = y)
        res <- rectGrob(x = coord$x,
            y = coord$y,
            width = coordX$size,
            height = unit(fontSize, "bigpts"),
            gp = gpar(fill = convertedAnnotations, col = borderColor))
    } else {
        a <- x
        x <- unit(1, "npc") - y
        y <- unit(1, "npc") - a

        coord <- expand.grid(y = y, x = x)
        res <- rectGrob(x = coord$x,
            y = coord$y,
            width = unit(fontSize, "bigpts"),
            height = coordX$size,
            gp = gpar(fill = convertedAnnotations, col = borderColor))
=======
draw_annotations = function(converted_annotations, border_color, gaps, fontsize, horizontal){
    n = ncol(converted_annotations)
    m = nrow(converted_annotations)

    coord_x = find_coordinates(m, gaps)

    x = coord_x$coord - 0.5 * coord_x$size

    # y = cumsum(rep(fontsize, n)) - 4 + cumsum(rep(2, n))
    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1
    y = unit(y, "bigpts")

    if(horizontal){
        coord = expand.grid(x = x, y = y)
        res = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = unit(fontsize, "bigpts"), gp = gpar(fill = converted_annotations, col = border_color))
    }
    else{
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a

        coord = expand.grid(y = y, x = x)
        res = rectGrob(x = coord$x, y = coord$y, width = unit(fontsize, "bigpts"), height = coord_x$size, gp = gpar(fill = converted_annotations, col = border_color))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
    }

    return(res)
}

<<<<<<< HEAD
.drawAnnotationNames <- function(annotations, fontSize, horizontal) {
    n <- ncol(annotations)

    x <- unit(3, "bigpts")

    y <- cumsum(rep(fontSize, n)) +
        cumsum(rep(2, n)) -
        fontSize / 2 + 1

    y <- unit(y, "bigpts")

    if (horizontal) {
        res <- textGrob(colnames(annotations),
            x = x,
            y = y,
            hjust = 0,
            gp = gpar(fontSize = fontSize, fontface = 2))
    } else {
        a <- x
        x <- unit(1, "npc") - y
        y <- unit(1, "npc") - a

        res <- textGrob(colnames(annotations),
            x = x,
            y = y,
            vjust = 0.5,
            hjust = 0,
            rot = 270,
            gp = gpar(fontSize = fontSize, fontface = 2))
=======
draw_annotation_names = function(annotations, fontsize, horizontal){
    n = ncol(annotations)

    x = unit(3, "bigpts")

    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1
    y = unit(y, "bigpts")

    if(horizontal){
        res = textGrob(colnames(annotations), x = x, y = y, hjust = 0, gp = gpar(fontsize = fontsize, fontface = 2))
    }
    else{
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a

        res = textGrob(colnames(annotations), x = x, y = y, vjust = 0.5, hjust = 0, rot = 270, gp = gpar(fontsize = fontsize, fontface = 2))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
    }

    return(res)
}

<<<<<<< HEAD
.drawAnnotationLegend <- function(annotation,
    annotationColors,
    borderColor,
    ...) {

    y <- unit(1, "npc")
    textHeight <- unit(1,
        "grobheight",
        textGrob("FGH", gp = gpar(...)))

    res <- gList()
    for (i in names(annotation)) {
        res[[i]] <- textGrob(i,
            x = 0,
            y = y,
            vjust = 1,
            hjust = 0,
            gp = gpar(fontface = "bold", ...))

        y <- y - 1.5 * textHeight
        if (is.character(annotation[[i]]) |
                is.factor(annotation[[i]])) {
            n <- length(annotationColors[[i]])
            yy <- y - (1:n - 1) * 2 * textHeight

            res[[paste(i, "r")]] <- rectGrob(x = unit(0, "npc"),
                y = yy,
                hjust = 0,
                vjust = 1,
                height = 2 * textHeight,
                width = 2 * textHeight,
                gp = gpar(col = borderColor, fill = annotationColors[[i]]))

            res[[paste(i, "t")]] <- textGrob(names(annotationColors[[i]]),
                x = textHeight * 2.4,
                y = yy - textHeight,
                hjust = 0,
                vjust = 0.5,
                gp = gpar(...))

            y <- y - n * 2 * textHeight
        } else {
            yy <- y - 8 * textHeight + seq(0, 1, 0.25)[-1] * 8 * textHeight
            h <- 8 * textHeight * 0.25

            res[[paste(i, "r")]] <- rectGrob(x = unit(0, "npc"),
                y = yy,
                hjust = 0,
                vjust = 1,
                height = h,
                width = 2 * textHeight,
                gp = gpar(col = NA,
                        fill = colorRampPalette(annotationColors[[i]])(4)))
            res[[paste(i, "r2")]] <- rectGrob(x = unit(0, "npc"),
                y = y,
                hjust = 0,
                vjust = 1,
                height = 8 * textHeight,
                width = 2 * textHeight,
                gp = gpar(col = borderColor, fill = NA))

            txt <- rev(range(grid.pretty(range(annotation[[i]],
                na.rm = TRUE))))

            yy <- y - c(1, 7) * textHeight
            res[[paste(i, "t")]] <- textGrob(txt,
                x = textHeight * 2.4,
                y = yy,
                hjust = 0,
                vjust = 0.5,
                gp = gpar(...))
            y <- y - 8 * textHeight
=======
draw_annotation_legend = function(annotation, annotation_colors, border_color, ...){
    y = unit(1, "npc")
    text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))

    res = gList()
    for(i in names(annotation)){
        res[[i]] = textGrob(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))

        y = y - 1.5 * text_height
        if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
            n = length(annotation_colors[[i]])
            yy = y - (1:n - 1) * 2 * text_height

            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = 2 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]]))
            res[[paste(i, "t")]] = textGrob(names(annotation_colors[[i]]), x = text_height * 2.4, y = yy - text_height, hjust = 0, vjust = 0.5, gp = gpar(...))

            y = y - n * 2 * text_height

        }
        else{
            yy = y - 8 * text_height + seq(0, 1, 0.25)[-1] * 8 * text_height
            h = 8 * text_height * 0.25

            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = 2 * text_height, gp = gpar(col = NA, fill = colorRampPalette(annotation_colors[[i]])(4)))
            res[[paste(i, "r2")]] = rectGrob(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = 8 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = NA))

            txt = rev(range(grid.pretty(range(annotation[[i]], na.rm = TRUE))))
            yy = y - c(1, 7) * text_height
            res[[paste(i, "t")]]  = textGrob(txt, x = text_height * 2.4, y = yy, hjust = 0, vjust = 0.5, gp = gpar(...))
            y = y - 8 * text_height
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
        }
        y <- y - 1.5 * textHeight
    }

<<<<<<< HEAD
    res <- gTree(children = res)
=======
    res = gTree(children = res)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

<<<<<<< HEAD
.drawMain <- function(text, ...) {
    res <- textGrob(text, gp = gpar(fontface = "bold", ...))
=======
draw_main = function(text, ...){
    res = textGrob(text, gp = gpar(fontface = "bold", ...))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    return(res)
}

vplayout <- function(x, y) {
    return(viewport(layout.pos.row = x, layout.pos.col = y))
}

<<<<<<< HEAD
.heatmapMotor <- function(matrix,
        borderColor,
        cellWidth,
        cellHeight,
        treeCol,
        treeRow,
        treeHeightCol,
        treeHeightRow,
        fileName,
        width,
        height,
        breaks,
        color,
        legend,
        annotationRow,
        annotationCol,
        annotationColors,
        annotationLegend,
        annotationNamesRow,
        annotationNamesCol,
        main,
        fontSize,
        fontSizeRow,
        fontSizeCol,
        fmat,
        fontSizeNumber,
        numberColor,
        gapsCol,
        gapsRow,
        labelsRow,
        labelsCol,
        ...) {
        # Set layout
        lo <- .lo(coln = labelsCol,
            rown = labelsRow,
            nrow = nrow(matrix),
            ncol = ncol(matrix),
            cellWidth = cellWidth,
            cellHeight = cellHeight,
            treeHeightCol = treeHeightCol,
            treeHeightRow = treeHeightRow,
            legend = legend,
            annotationCol = annotationCol,
            annotationRow = annotationRow,
            annotationColors = annotationColors,
            annotationLegend = annotationLegend,
            annotationNamesRow = annotationNamesRow,
            annotationNamesCol = annotationNamesCol,
            main = main,
            fontSize = fontSize,
            fontSizeRow = fontSizeRow,
            fontSizeCol = fontSizeCol,
            gapsRow = gapsRow,
            gapsCol = gapsCol,
            ...)

        res <- lo$gt
        mindim <- lo$mindim

        if (!is.na(fileName)) {
            if (is.na(height)) {
                height <- convertHeight(gtable_height(res),
                    "inches",
                    valueOnly = TRUE)
            }
            if (is.na(width)) {
                width <- convertWidth(gtable_width(res),
                    "inches",
                    valueOnly = TRUE)
            }

            # Get file type
            r <- regexpr("\\.[a-zA-Z]*$", fileName)
            if (r == -1)
                stop("Improper fileName")
            ending <- substr(fileName,
                r + 1,
                r + attr(r, "match.length"))

            f <- switch(ending,
                pdf = function(x, ...)
                pdf(x, ...),
                png = function(x, ...)
                png(x, units = "in",
                    res = 300, ...),
                jpeg = function(x, ...)
                jpeg(x, units = "in",
                    res = 300, ...),
                jpg = function(x, ...)
                jpeg(x, units = "in",
                    res = 300, ...),
                tiff = function(x, ...)
                tiff(x,
                    units = "in",
                    res = 300,
                    compression = "lzw",
                    ...),
                bmp = function(x, ...)
                bmp(x, units = "in",
                    res = 300, ...),
                stop("File type should be: pdf, png, bmp, jpg, tiff"))

            # print(sprintf("height:%f width:%f", height, width))

            # gt = .heatmapMotor(matrix,
            #     cellWidth = cellWidth,
            #     cellHeight = cellHeight,
            #     borderColor = borderColor,
            #     treeCol = treeCol,
            #     treeRow = treeRow,
            #     treeHeightCol = treeHeightCol,
            #     treeHeightRow = treeHeightRow,
            #     breaks = breaks,
            #     color = color,
            #     legend = legend,
            #     annotationCol = annotationCol,
            #     annotationRow = annotationRow,
            #     annotationColors = annotationColors,
            #     annotationLegend = annotationLegend,
            #     fileName = NA, main = main,
            #     fontSize = fontSize,
            #     fontSizeRow = fontSizeRow,
            #     fontSizeCol = fontSizeCol,
            #     fmat = fmat,
            #     fontSizeNumber =  fontSizeNumber,
            #     numberColor = numberColor,
            #     labelsRow = labelsRow,
            #     labelsCol = labelsCol,
            #     gapsCol = gapsCol,
            #     gapsRow = gapsRow, ...)

            f(fileName, height = height, width = width)
            gt <- .heatmapMotor(matrix,
                cellWidth = cellWidth,
                cellHeight = cellHeight,
                borderColor = borderColor,
                treeCol = treeCol,
                treeRow = treeRow,
                treeHeightCol = treeHeightCol,
                treeHeightRow = treeHeightRow,
                breaks = breaks,
                color = color,
                legend = legend,
                annotationCol = annotationCol,
                annotationRow = annotationRow,
                annotationColors = annotationColors,
                annotationLegend = annotationLegend,
                annotationNamesRow = annotationNamesRow,
                annotationNamesCol = annotationNamesCol,
                fileName = NA,
                main = main,
                fontSize = fontSize,
                fontSizeRow = fontSizeRow,
                fontSizeCol = fontSizeCol,
                fmat = fmat,
                fontSizeNumber = fontSizeNumber,
                numberColor = numberColor,
                labelsRow = labelsRow,
                labelsCol = labelsCol,
                gapsCol = gapsCol,
                gapsRow = gapsRow,
                ...)
            grid.draw(gt)
            dev.off()

            return(gt)
=======
heatmap_motor = function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, annotation_names_row, annotation_names_col, main, fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, number_color, gaps_col, gaps_row, labels_row, labels_col, ...){
    # Set layout
    lo = lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, gaps_row = gaps_row, gaps_col = gaps_col,  ...)

    res = lo$gt
    mindim = lo$mindim

    if(!is.na(filename)){
        if(is.na(height)){
            height = convertHeight(gtable_height(res), "inches", valueOnly = TRUE)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
        }

        # Omit border color if cell size is too small
        if (mindim < 3){
            borderColor <- NA
        }

<<<<<<< HEAD
        # Draw title
        if (!is.na(main)) {
            elem <- .drawMain(main, fontSize = 1.3 * fontSize, ...)
            res <- gtable_add_grob(res,
                elem,
                t = 1,
                l = 3,
                name = "main",
                clip = "off")
        }

        # Draw tree for the columns
        if (!.is.na2(treeCol) & treeHeightCol != 0) {
            elem <- .drawDendrogram(treeCol, gapsCol, horizontal = TRUE)
            res <- gtable_add_grob(res,
                elem,
                t = 2,
                l = 3,
                name = "col_tree")
        }

        # Draw tree for the rows
        if (!.is.na2(treeRow) & treeHeightRow != 0) {
            elem <- .drawDendrogram(treeRow, gapsRow, horizontal = FALSE)
            res <- gtable_add_grob(res,
                elem,
                t = 4,
                l = 1,
                name = "row_tree")
        }

        # Draw matrix
        elem <- .drawMatrix(matrix,
            borderColor,
            gapsRow,
            gapsCol,
            fmat,
            fontSizeNumber,
            numberColor)

        res <- gtable_add_grob(res,
                elem,
                t = 4,
                l = 3,
                clip = "off",
                name = "matrix")

        # Draw colnames
        if (length(labelsCol) != 0) {
            pars <- list(labelsCol,
                gaps = gapsCol,
                fontSize = fontSizeCol,
                ...)
            elem <- do.call(.drawColnames, pars)
            res <- gtable_add_grob(res,
                elem,
                t = 5,
                l = 3,
                clip = "off",
                name = "col_names")
        }

        # Draw rownames
        if (length(labelsRow) != 0) {
            pars <- list(labelsRow,
                gaps = gapsRow,
                fontSize = fontSizeRow, ...)
            elem <- do.call(.drawRownames, pars)
            res <- gtable_add_grob(res,
                elem,
                t = 4,
                l = 4,
                clip = "off",
                name = "row_names")
        }

        # Draw annotation tracks on cols
        if (!.is.na2(annotationCol)) {
            # Draw tracks
            convertedAnnotation <- .convertAnnotations(annotationCol,
                annotationColors)
            elem <- .drawAnnotations(convertedAnnotation,
                borderColor,
                gapsCol,
                fontSize,
                horizontal = TRUE)
            res <- gtable_add_grob(res,
                elem,
                t = 3,
                l = 3,
                clip = "off",
                name = "col_annotation")

            # Draw names
            if (annotationNamesCol) {
                elem <- .drawAnnotationNames(annotationCol,
                    fontSize,
                    horizontal = TRUE)
                res <- gtable_add_grob(res,
                    elem,
                    t = 3,
                    l = 4,
                    clip = "off",
                    name = "col_annotation_names")
            }
        }

        # Draw annotation tracks on rows
        if (!.is.na2(annotationRow)) {
            # Draw tracks
            convertedAnnotation <- .convertAnnotations(annotationRow,
                annotationColors)
            elem <- .drawAnnotations(convertedAnnotation,
                borderColor,
                gapsRow,
                fontSize,
                horizontal = FALSE)
            res <- gtable_add_grob(res,
                elem,
                t = 4,
                l = 2,
                clip = "off",
                name = "row_annotation")

            # Draw names
            if (annotationNamesRow) {
                elem <- .drawAnnotationNames(annotationRow,
                    fontSize,
                    horizontal = FALSE)
                res <- gtable_add_grob(res,
                    elem,
                    t = 5,
                    l = 2,
                    clip = "off",
                    name = "row_annotation_names")
            }
        }

        # Draw annotation legend
        annotation <- c(annotationCol[length(annotationCol):1],
            annotationRow[length(annotationRow):1])
        annotation <- annotation[unlist(lapply(annotation,
            function(x) !.is.na2(x)))]

        if (length(annotation) > 0 & annotationLegend) {
            elem <- .drawAnnotationLegend(annotation,
                annotationColors,
                borderColor,
                fontSize = fontSize,
                ...)

            t <- ifelse(is.null(labelsRow), 4, 3)
            res <- gtable_add_grob(res,
                elem,
                t = t,
                l = 6,
                b = 5,
                clip = "off",
                name = "annotationLegend")
        }
=======
        # Get file type
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if(r == -1) stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))

        f = switch(ending,
            pdf = function(x, ...) pdf(x, ...),
            png = function(x, ...) png(x, units = "in", res = 300, ...),
            jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
            jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
            tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
            bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
            stop("File type should be: pdf, png, bmp, jpg, tiff")
        )

        # print(sprintf("height:%f width:%f", height, width))

        # gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)

        f(filename, height = height, width = width)
        gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)
        grid.draw(gt)
        dev.off()

        return(gt)
    }

    # Omit border color if cell size is too small
    if(mindim < 3) border_color = NA

    # Draw title
    if(!is.na(main)){
        elem = draw_main(main, fontsize = 1.3 * fontsize, ...)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", clip = "off")
    }

    # Draw tree for the columns
    if(!is.na2(tree_col) & treeheight_col != 0){
        elem = draw_dendrogram(tree_col, gaps_col, horizontal = TRUE)
        res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }

    # Draw tree for the rows
    if(!is.na2(tree_row) & treeheight_row != 0){
        elem = draw_dendrogram(tree_row, gaps_row, horizontal = FALSE)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }

    # Draw matrix
    elem = draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", name = "matrix")

    # Draw colnames
    if(length(labels_col) != 0){
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, ...)
        elem = do.call(draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "col_names")
    }

    # Draw rownames
    if(length(labels_row) != 0){
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, ...)
        elem = do.call(draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 4, clip = "off", name = "row_names")
    }

    # Draw annotation tracks on cols
    if(!is.na2(annotation_col)){
        # Draw tracks
        converted_annotation = convert_annotations(annotation_col, annotation_colors)
        elem = draw_annotations(converted_annotation, border_color, gaps_col, fontsize, horizontal = TRUE)
        res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", name = "col_annotation")

        # Draw names
        if(annotation_names_col){
            elem = draw_annotation_names(annotation_col, fontsize, horizontal = TRUE)
            res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", name = "col_annotation_names")
        }
    }

    # Draw annotation tracks on rows
    if(!is.na2(annotation_row)){
        # Draw tracks
        converted_annotation = convert_annotations(annotation_row, annotation_colors)
        elem = draw_annotations(converted_annotation, border_color, gaps_row, fontsize, horizontal = FALSE)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", name = "row_annotation")

        # Draw names
        if(annotation_names_row){
            elem = draw_annotation_names(annotation_row, fontsize, horizontal = FALSE)
            res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", name = "row_annotation_names")
        }
    }

    # Draw annotation legend
    annotation = c(annotation_col[length(annotation_col):1], annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]

    if(length(annotation) > 0 & annotation_legend){
        elem = draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)

        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, clip = "off", name = "annotation_legend")
    }

    # Draw legend
    if(!is.na2(legend)){
        elem = draw_legend(color, breaks, legend, fontsize = fontsize, ...)

        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, clip = "off", name = "legend")
    }

    return(res)
}
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

        # Draw legend
        if (!.is.na2(legend)) {
            elem <- .drawLegend(color, breaks, legend, fontSize = fontSize, ...)

            t <- ifelse(is.null(labelsRow), 4, 3)
            res <- gtable_add_grob(res,
                    elem,
                    t = t,
                    l = 5,
                    b = 5,
                    clip = "off",
                    name = "legend")
        }

        return(res)
    }

.generateBreaks <- function(x, n, center = FALSE) {
    if (center) {
        m <- max(abs(c(min(x,na.rm = TRUE),
            max(x, na.rm = TRUE))))
        res <- seq(-m, m, length.out = n + 1)
    } else {
        res <- seq(min(x, na.rm = TRUE),
            max(x, na.rm = TRUE),
            length.out = n + 1)
    }

    return(res)
}

.scaleVecColours <- function(x, col = rainbow(10), breaks = NA) {
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = TRUE))])
}

.scaleColours <- function(mat,
    col = rainbow(10),
    breaks = NA) {
    mat <- as.matrix(mat)
    return(matrix(
        .scaleVecColours(as.vector(mat), col = col, breaks = breaks),
            nrow(mat),
            ncol(mat),
            dimnames = list(rownames(mat), colnames(mat))))
}

<<<<<<< HEAD
## changed the original clusterMat() in the pheatmap.r
.clusterMat <- function(mat, labels, distance, method) {
    # this funciton is going to change the .clusterMat() in pheatmap

    if (!(method %in% c("ward.D",
            "ward.D2",
            "ward",
            "single",
            "complete",
            "average",
            "mcquitty",
            "median",
            "centroid"))) {
        stop("clustering method has to one form the list:
            'ward',
            'ward.D',
            'ward.D2',
            'single',
            'complete',
            'average',
            'mcquitty',
            'median'
            or 'centroid'.")
    }

    class.label <- unique(labels)

    nGroup <- length(class.label) # [#group]
    # get "hclust" object for each group then wrap them up as group.hclust

    # distance function preparation
    dis <- function(mat, distance) {
        if (!(distance[1] %in% c("correlation",
                "euclidean",
                "maximum",
                "manhattan",
                "canberra",
                "binary",
                "minkowski")) &
                !methods::is(distance, "dist")) {
            stop("distance has to be a dissimilarity structure as produced by
                dist or one measure  form the list: 'correlation', 'euclidean',
                'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
        }

        if (distance[1] == "correlation") {
            # this part should be confirmed whether being wrong?
            #ToDo: how is the correlation matrix converted to a dsit matrix
            d <- stats::as.dist(1 - stats::cor(t(mat)))
        } else {
            d <- stats::dist(mat, method = distance)
        }

        return(d)
    }

    # initiate the final returning value: a "hclust" object
    cum.hclust <- list()

    if (nGroup == 1) {
        # matrix has only 1 group
        if (length(labels) == 1) {
            stop("only one row/column for the matrix")
        }
        group.hclust <- stats::hclust(dis(mat = mat,
                distance = distance),
            method = method)

        cum.hclust <- group.hclust
    } else {
        #  matrix has more than 1 groups
        group.hclust <- sapply(class.label, function(x) {
            # get the positions of class label
            class.pos <- which(labels == x)

            if (length(class.pos) == 1) {
                # if only 1 row in the group return a manually made "hclust" object
                sub.hclust <- as.list(1:7)
                names(sub.hclust) <- c("merge",
                    "height",
                    "order",
                    "labels",
                    "method",
                    "call",
                    "dist.method")

                class(sub.hclust) <- "hclust"
                sub.hclust$merge <- matrix(c(0, 0), nrow = 1)
                sub.hclust$height <- 0
                sub.hclust$order <- 1
                return(sub.hclust)
            } else if (length(class.pos) > 1) {
                # if >1 rows return the "hclust" object
                return(stats::hclust(dis(mat = mat[class.pos,],
                    distance = distance),
                    method = method))
            }
        })
        # the length(group.hclust) is the [#group] == nGroup   ,
        # group.hclust[[i]] to get each "hclust"

        # then modify the "hclust" object and get them merged into one
        # "hclust" object

        # initiate the final "hclust" object
        cum.hclust <- group.hclust[, nGroup]

        # merge function preparation
        mergeHclust <- function(hclust1, hclust2) {
                # "hclust" object modifying function preparation
                if (hclust1$merge[1, 1] == 0 &
                        hclust2$merge[1, 1] == 0) {
                    # both groups have only 1 row
                    hclustCom <- as.list(seq(1,7))
                    names(hclustCom) <-
                        c("merge",
                        "height",
                        "order",
                        "labels",
                        "method",
                        "call",
                        "dist.method")

                    class(hclustCom) <- "hclust"
                    hclustCom$merge <- matrix(c(-1,-2), nrow = 1)
                    # check for different matrix whether 1 should be good
                    hclustCom$height <- 1
                    hclustCom$order <- c(1, 2)
                    return(hclustCom)
                } else if (hclust1$merge[1, 1] != 0 &
                        hclust2$merge[1, 1] != 0) {
                    # both group have >1 rows

                    # nodes in the hclust1 group, so actually the #rows should
                    # be dim()[1]+1
                    row.1 <- dim(hclust1$merge)[1]
                    # nodes in the hclust2 group
                    row.2 <- dim(hclust2$merge)[1]
                    hclustCom <- list()
                    mer <- hclust2$merge
                    # modify the hclust2$merge matrix
                    hclustCom$merge <- (mer > 0) *
                        (mer + row.1) + (mer < 0) *
                        (mer - row.1 - 1)
                    # combine the merge matrix from the 2 groups
                    hclustCom$merge <- rbind(hclust1$merge,
                        hclustCom$merge)
                    hclustCom$merge <- rbind(hclustCom$merge,
                        c(row.1, row.1 + row.2))
                    hclustCom$height <- c(hclust1$height, hclust2$height)
                    # check for different matrix whether 1 should be good
                    hclustCom$height <- c(hclustCom$height,
                        max(hclustCom$height) + 1)
                    hclustCom$order <- c(hclust1$order,
                        hclust2$order + row.1 + 1)
                    class(hclustCom) <- "hclust"
                } else {
                    # one group has only 1 row, the other group has >1 rows
                    if (hclust1$merge[1, 1] == 0) {
                        # hclust1 has 1 row , hclust2 has >1 rows

                        # nodes in the hclust2 group
                        row.2 <- dim(hclust2$merge)[1]
                        hclustCom <- as.list(1:7)
                        names(hclustCom) <- c("merge",
                                "height",
                                "order",
                                "labels",
                                "method",
                                "call",
                                "dist.method")
                        class(hclustCom) <- "hclust"
                        mer <- hclust2$merge
                        hclustCom$merge <- (mer > 0) *
                            (mer) +
                            (mer < 0) *
                            (mer - 1)
                        hclustCom$merge <- rbind(hclustCom$merge,
                            c(-1, row.2))
                        # check for different matrix whether 1 should be good
                        hclustCom$height <- c(hclust2$height,
                            max(hclust2$height) + 1)
                        hclustCom$order <- c(1, hclust2$order + 1)
                    } else if (hclust2$merge[1, 1] == 0) {
                        # the hclust1 has >1 rows , and hclust2 has 1 row

                        #nodes in the hclust1 group
                        row.1 <- dim(hclust1$merge)[1]
                        hclustCom <- as.list(seq(1,7))
                        names(hclustCom) <-
                            c("merge",
                                "height",
                                "order",
                                "labels",
                                "method",
                                "call",
                                "dist.method")
                        class(hclustCom) <- "hclust"
                        hclustCom$merge <- hclust1$merge
                        hclustCom$merge <- rbind(hclustCom$merge,
                            c(row.1,- (row.1 + 2)))
                        hclustCom$height <- c(hclust1$height,
                            max(hclust1$height) + 1)
                        hclustCom$order <- c(hclust1$order,
                            max(hclust1$order) + 1)
                    }
                }

                return(hclustCom)
            }

        # merge the "hclust" object into the final one "hclust" object
        for (i in seq(nGroup - 1, 1, -1)) {
            cum.hclust <- mergeHclust(group.hclust[, i], cum.hclust)
        }
    }

    cum.hclust$labels <- NULL
    cum.hclust$call <- NULL
    cum.hclust$method <- NULL
    cum.hclust$dist.method <- NULL

    return(cum.hclust)
=======
## changed the original cluster_mat() in the pheatmap.r
cluster_mat <- function(mat, labels, distance, method){
  # this funciton is going to change the cluster_mat() in pheatmap

  if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
    stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }

  class.label <- unique(labels)

  nGroup <- length(class.label)    # [#group]
  # get "hclust" object for each group then wrap them up as group.hclust

  # distance function preparation
  dis <- function(mat, distance) {
    if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & !methods::is(distance, "dist")){
      stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if(distance[1] == "correlation"){
      # this part should be confirmed whether being wrong? ToDo: how is the correlation matrix converted to a dsit matrix
      d <- stats::as.dist(1 - stats::cor(t(mat)))
    }
    else{
      d <- stats::dist(mat, method = distance)
    }

    return(d)
  }

  # initiate the final returning value: a "hclust" object
  cum.hclust <- list()

  if(nGroup==1){  # matrix has only 1 group
    if(length(labels)==1){
      stop("only one row/column for the matrix")
    }
    group.hclust <-  stats::hclust(dis(mat = mat, distance = distance), method = method)

    cum.hclust <- group.hclust

  }else {           #  matrix has more than 1 groups
    group.hclust <- sapply(class.label, function(x) {

      class.pos <- which(labels==x)   # get the positions of class label
      if(length(class.pos)==1){  # if only 1 row in the group return a manually made "hclust" object
        sub.hclust <- as.list(1:7)
        names(sub.hclust) <- c("merge","height", "order", "labels","method","call","dist.method")
        class(sub.hclust) <- "hclust"
        sub.hclust$merge <- matrix(c(0,0), nrow = 1)
        sub.hclust$height <- 0
        sub.hclust$order <- 1
        return(sub.hclust)
      }else if(length(class.pos)>1){   # if >1 rows return the "hclust" object
        return(stats::hclust(dis(mat =  mat[class.pos,], distance = distance ),method = method ))
      } }
    )
    # the length(group.hclust) is the [#group] == nGroup   ,   group.hclust[[i]] to get each "hclust"

    # then modify the "hclust" object and get them merged into one "hclust" object
    cum.hclust <- group.hclust[,nGroup]   # initiate the final "hclust" object

    # merge function preparation
    merge_hclust <- function(hclust1, hclust2){   # "hclust" object modifying function preparation
      if(hclust1$merge[1,1]==0 & hclust2$merge[1,1]==0){ # both groups have only 1 row
        hclust.com <- as.list(1:7)
        names(hclust.com) <- c("merge","height", "order", "labels","method","call","dist.method")
        class(hclust.com) <- "hclust"
        hclust.com$merge <- matrix(c(-1,-2), nrow = 1)
        hclust.com$height <- 1   # check for different matrix whether 1 should be good
        hclust.com$order <- c(1,2)
        return(hclust.com)
      }else if( hclust1$merge[1,1]!=0 & hclust2$merge[1,1]!=0) {  # both group have >1 rows
        row.1 <- dim(hclust1$merge)[1]  #   #nodes in the hclust1 group, so actually the #rows shouls be dim()[1]+1
        row.2 <- dim(hclust2$merge)[1]  #   #nodes in the hclust2 group
        hclust.com <- list()
        mer <- hclust2$merge
        hclust.com$merge <- (mer>0) * (mer+row.1) + (mer<0) * (mer-row.1-1)  # modify the hclust2$merge matrix
        hclust.com$merge <- rbind(hclust1$merge, hclust.com$merge)     # combine the merge matrix from the 2 groups
        hclust.com$merge <- rbind(hclust.com$merge, c(row.1, row.1+row.2))
        hclust.com$height <- c(hclust1$height, hclust2$height)
        hclust.com$height <- c(hclust.com$height, max(hclust.com$height)+ 1 )  # check for different matrix whether 1 should be good
        hclust.com$order <- c(hclust1$order, hclust2$order+row.1+1)
        class(hclust.com) <- "hclust"
      }else{  # one group has only 1 row, the other group has >1 rows
        if(hclust1$merge[1,1]==0){   #  hclust1 has 1 row , hclust2 has >1 rows
          row.2 <- dim(hclust2$merge)[1]  #   #nodes in the hclust2 group
          hclust.com <- as.list(1:7)
          names(hclust.com) <- c("merge","height", "order", "labels","method","call","dist.method")
          class(hclust.com) <- "hclust"
          mer <- hclust2$merge
          hclust.com$merge <- (mer>0) * (mer) + (mer<0) * (mer-1)
          hclust.com$merge <- rbind(hclust.com$merge, c(-1, row.2))
          hclust.com$height <- c(hclust2$height, max(hclust2$height)+1)  # check for different matrix whether 1 should be good
          hclust.com$order <- c(1, hclust2$order+1)
        }else if(hclust2$merge[1,1]==0){   # the hclust1 has >1 rows , and hclust2 has 1 row
          row.1 <- dim(hclust1$merge)[1]  #   #nodes in the hclust1 group
          hclust.com <- as.list(1:7)
          names(hclust.com) <- c("merge","height", "order", "labels","method","call","dist.method")
          class(hclust.com) <- "hclust"
          hclust.com$merge <- hclust1$merge
          hclust.com$merge <- rbind(hclust.com$merge, c(row.1, -(row.1+2)))
          hclust.com$height <- c(hclust1$height, max(hclust1$height)+1)
          hclust.com$order <- c(hclust1$order, max(hclust1$order)+1)
        }

      }

      return(hclust.com)
    }

    # merge the "hclust" object into the final one "hclust" object
    for(i in seq(nGroup-1, 1,-1)){
      cum.hclust <- merge_hclust( group.hclust[,i], cum.hclust)
    }

  }

  cum.hclust$labels <-NULL
  cum.hclust$call <- NULL
  cum.hclust$method <- NULL
  cum.hclust$dist.method <- NULL

  return(cum.hclust)

>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
}

.scaleRows <- function(x) {
    m <- base::apply(x, 1, mean, na.rm = TRUE)
    s <- base::apply(x, 1, stats::sd, na.rm = TRUE)
    return((x - m) / s)
}

.scaleMat <- function(mat, scale) {
    if (!(scale %in% c("none", "row", "column"))) {
        stop("scale argument shoud take values: 'none', 'row' or 'column'")
    }
    mat <- switch(scale,
            none = mat,
            row = .scaleRows(mat),
            column = t(.scaleRows(t(mat))))
    return(mat)
}

.generateAnnotationColours <- function(annotation,
    annotationColors,
    drop) {

        if (.is.na2(annotationColors)) {
            annotationColors <- list()
        }
<<<<<<< HEAD
        count <- 0
        for (i in 1:length(annotation)) {
            annotation[[i]] <- annotation[[i]][!is.na(annotation[[i]])]
            if (is.character(annotation[[i]]) |
                    is.factor(annotation[[i]])) {
                if (is.factor(annotation[[i]]) & !drop) {
                    count <- count + length(levels(annotation[[i]]))
                } else {
                    count <- count + length(unique(annotation[[i]]))
                }
=======
    }

    factor_colors = dscale(factor(1:count), hue_pal(l = 75))

    cont_counter = 2
    for(i in 1:length(annotation)){
        if(!(names(annotation)[i] %in% names(annotation_colors))){
            if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
                n = length(unique(annotation[[i]]))
                if (is.factor(annotation[[i]]) & !drop){
                    n = length(levels(annotation[[i]]))
                }
                ind = sample(1:length(factor_colors), n)
                annotation_colors[[names(annotation)[i]]] = factor_colors[ind]
                l = levels(as.factor(annotation[[i]]))
                l = l[l %in% unique(annotation[[i]])]
                if (is.factor(annotation[[i]]) & !drop){
                    l = levels(annotation[[i]])
                }

                names(annotation_colors[[names(annotation)[i]]]) = l
                factor_colors = factor_colors[-ind]
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
            }
        }

        factorColors <- scales::dscale(factor(seq(1,count)), hue_pal(l = 75))

        contCounter <- 2
        for (i in 1:length(annotation)) {
            if (!(names(annotation)[i] %in% names(annotationColors))) {
                if (is.character(annotation[[i]]) |
                        is.factor(annotation[[i]])) {
                    n <- length(unique(annotation[[i]]))

                    if (is.factor(annotation[[i]]) & !drop) {
                        n <- length(levels(annotation[[i]]))
                    }

                    ind <- sample(1:length(factorColors), n)
                    annotationColors[[names(annotation)[i]]] <-
                        factorColors[ind]
                    l <- levels(as.factor(annotation[[i]]))
                    l <- l[l %in% unique(annotation[[i]])]
                    if (is.factor(annotation[[i]]) & !drop) {
                        l <- levels(annotation[[i]])
                    }

                    names(annotationColors[[names(annotation)[i]]]) <- l
                    factorColors <- factorColors[-ind]
                } else {
                    annotationColors[[names(annotation)[i]]] <-
                        brewer_pal("seq", contCounter)(5)[1:4]
                    contCounter <- contCounter + 1
                }
            }
        }
        return(annotationColors)
    }

<<<<<<< HEAD

.findGaps <- function(tree, cutreeN) {
    v <- stats::cutree(tree, cutreeN)[tree$order]
    gaps <- which((v[-1] - v[-length(v)]) != 0)
=======
find_gaps = function(tree, cutree_n){
    v = stats::cutree(tree, cutree_n)[tree$order]
    gaps = which((v[-1] - v[-length(v)]) != 0)

>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
}

.is.na2 <- function(x) {
    if (is.list(x) | length(x) > 1) {
        return(FALSE)
    }
    if (length(x) == 0) {
        return(TRUE)
    }

    return(is.na(x))
}

.identity2 <- function(x, ...) {
    return(x)
}

<<<<<<< HEAD
#' @title A function to draw clustered heatmaps.
#' @description A function to draw clustered heatmaps where one has better
#'  control over some graphical parameters such as cell size, etc.
=======
#' A function to draw clustered heatmaps.
#'
#' A function to draw clustered heatmaps where one has better control over some graphical
#' parameters such as cell size, etc.
#'
#' The function also allows to aggregate the rows using kmeans clustering. This is
#' advisable if number of rows is so big that R cannot handle their hierarchical
#' clustering anymore, roughly more than 1000. Instead of showing all the rows
#' separately one can cluster the rows in advance and show only the cluster centers.
#' The number of clusters can be tuned with parameter kmeans_k.
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
#'
#' The function also allows to aggregate the rows using kmeans clustering.
#'  This is advisable if number of rows is so big that R cannot handle their
#'  hierarchical clustering anymore, roughly more than 1000. Instead of showing
#'  all the rows separately one can cluster the rows in advance and show only
#'  the cluster centers. The number of clusters can be tuned with parameter
#'  kmeansK.
#' @param mat numeric matrix of the values to be plotted.
#' @param color vector of colors used in heatmap.
<<<<<<< HEAD
#' @param kmeansK the number of kmeans clusters to make, if we want to
#'  agggregate the rows before drawing heatmap. If NA then the rows are not
#'  aggregated.
#' @param breaks Numeric vector. A sequence of numbers that covers the range
#'  of values in the normalized `counts`. Values in the normalized `matrix` are
#'  assigned to each bin in `breaks`. Each break is assigned to a unique color
#'  from `col`. If NULL, then breaks are calculated automatically. Default NULL.
#' @param borderColor color of cell borders on heatmap, use NA if no border
#'  should be drawn.
#' @param cellWidth individual cell width in points. If left as NA, then the
#'  values depend on the size of plotting window.
#' @param cellHeight individual cell height in points. If left as NA, then the
#'  values depend on the size of plotting window.
#' @param scale character indicating if the values should be centered and
#'  scaled in either the row direction or the column direction, or none.
#'  Corresponding values are \code{"row"}, \code{"column"} and \code{"none"}.
#' @param clusterRows boolean values determining if rows should be clustered or
#'  \code{hclust} object,
#' @param clusterCols boolean values determining if columns should be clustered
#'  or \code{hclust} object.
#' @param clusteringDistanceRows distance measure used in clustering rows.
#'  Possible values are \code{"correlation"} for Pearson correlation and all
#'  the distances supported by \code{\link{dist}}, such as \code{"euclidean"},
#'  etc. If the value is none of the above it is assumed that a distance matrix
#'  is provided.
#' @param clusteringDistanceCols distance measure used in clustering columns.
#'  Possible values the same as for clusteringDistanceRows.
#' @param clusteringMethod clustering method used. Accepts the same values as
#'  \code{\link{hclust}}.
#' @param clusteringCallback callback function to modify the clustering. Is
#'  called with two parameters: original \code{hclust} object and the matrix
#'  used for clustering. Must return a \code{hclust} object.
#' @param cutreeRows number of clusters the rows are divided into, based on the
#'  hierarchical clustering (using cutree), if rows are not clustered, the
#'  argument is ignored
#' @param cutreeCols similar to \code{cutreeRows}, but for columns
#' @param treeHeightRow the height of a tree for rows, if these are clustered.
#'  Default value 50 points.
#' @param treeHeightCol the height of a tree for columns, if these are
#'  clustered. Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legendBreaks vector of breakpoints for the legend.
#' @param legendLabels vector of labels for the \code{legendBreaks}.
#' @param annotationRow data frame that specifies the annotations shown on left
#'  side of the heatmap. Each row defines the features for a specific row. The
#'  rows in the data and in the annotation are matched using corresponding row
#'  names. Note that color schemes takes into account if variable is continuous
#'  or discrete.
#' @param annotationCol similar to annotationRow, but for columns.
#' @param annotation deprecated parameter that currently sets the annotationCol
#'  if it is missing.
#' @param annotationColors list for specifying annotationRow and
#'  annotationCol track colors manually. It is  possible to define the colors
#'  for only some of the features. Check examples for  details.
#' @param annotationLegend boolean value showing if the legend for annotation
#'  tracks should be drawn.
#' @param annotationNamesRow boolean value showing if the names for row
#'  annotation tracks should be drawn.
#' @param annotationNamesCol boolean value showing if the names for column
#'  annotation tracks should be drawn.
#' @param dropLevels logical to determine if unused levels are also shown in
#'  the legend.
#' @param showRownames boolean specifying if column names are be shown.
#' @param showColnames boolean specifying if column names are be shown.
#' @param main the title of the plot
#' @param fontSize base fontsize for the plot
#' @param fontSizeRow fontsize for rownames (Default: fontsize)
#' @param fontSizeCol fontsize for colnames (Default: fontsize)
#' @param displayNumbers logical determining if the numeric values are also
#'  printed to the cells. If this is a matrix (with same dimensions as original
#'  matrix), the contents of the matrix are shown instead of original values.
#' @param numberFormat format strings (C printf style) of the numbers shown in
#'  cells. For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}"
#'  shows exponential notation (see more in \code{\link{sprintf}}).
#' @param numberColor color of the text
#' @param fontSizeNumber fontsize of the numbers displayed in cells
#' @param gapsRow vector of row indices that show shere to put gaps into
#'  heatmap. Used only if the rows are not clustered. See \code{cutreeRow}
#'  to see how to introduce gaps to clustered rows.
#' @param gapsCol similar to gapsRow, but for columns.
#' @param labelsRow custom labels for rows that are used instead of rownames.
#' @param labelsCol similar to labelsRow, but for columns.
#' @param fileName file path where to save the picture. Filetype is decided by
#'  the extension in the path. Currently following formats are supported: png,
#'  pdf, tiff, bmp, jpeg. Even if the plot does not fit into the plotting
#'  window, the file size is calculated so that the plot would fit there,
#'  unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param height manual option for determining the output file height in inches.
#' @param silent do not draw the plot (useful when using the gtable output)
#' @param rowLabel row cluster labels for semi-clustering
#' @param colLabel column cluster labels for semi-clustering
#' @param \dots graphical parameters for the text used in plot. Parameters
#'  passed to \code{\link{grid.text}}, see \code{\link{gpar}}.
#' @return
#' Invisibly a list of components
#' \itemize{
#'     \item \code{treeRow} the clustering of rows as \code{\link{hclust}}
#'       object
#'     \item \code{treeCol} the clustering of columns as \code{\link{hclust}}
#'       object
#'     \item \code{kmeans} the kmeans clustering of rows if parameter
#'       \code{kmeansK} was specified
#' }
=======
#' @param kmeans_k the number of kmeans clusters to make, if we want to agggregate the
#' rows before drawing heatmap. If NA then the rows are not aggregated.
#' @param breaks Numeric vector. A sequence of numbers that covers the range of values in the normalized `counts`. Values in the normalized `matrix` are assigned to each bin in `breaks`. Each break is assigned to a unique color from `col`. If NULL, then breaks are calculated automatically. Default NULL.
#' @param border_color color of cell borders on heatmap, use NA if no border should be
#' drawn.
#' @param cellwidth individual cell width in points. If left as NA, then the values
#' depend on the size of plotting window.
#' @param cellheight individual cell height in points. If left as NA,
#' then the values depend on the size of plotting window.
#' @param scale character indicating if the values should be centered and scaled in
#' either the row direction or the column direction, or none. Corresponding values are
#' \code{"row"}, \code{"column"} and \code{"none"}
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object,
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible
#' values are \code{"correlation"} for Pearson correlation and all the distances
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none
#' of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering columns. Possible
#' values the same as for clustering_distance_rows.
#' @param clustering_method clustering method used. Accepts the same values as
#' \code{\link{hclust}}.
#' @param clustering_callback callback function to modify the clustering. Is
#' called with two parameters: original \code{hclust} object and the matrix
#' used for clustering. Must return a \code{hclust} object.
#' @param cutree_rows number of clusters the rows are divided into, based on the
#'  hierarchical clustering (using cutree), if rows are not clustered, the
#' argument is ignored
#' @param cutree_cols similar to \code{cutree_rows}, but for columns
#' @param treeheight_row the height of a tree for rows, if these are clustered.
#' Default value 50 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered.
#' Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legend_breaks vector of breakpoints for the legend.
#' @param legend_labels vector of labels for the \code{legend_breaks}.
#' @param annotation_row data frame that specifies the annotations shown on left
#'  side of the heatmap. Each row defines the features for a specific row. The
#' rows in the data and in the annotation are matched using corresponding row
#'  names. Note that color schemes takes into account if variable is continuous
#'  or discrete.
#' @param annotation_col similar to annotation_row, but for columns.
#' @param annotation deprecated parameter that currently sets the annotation_col if it is missing
#' @param annotation_colors list for specifying annotation_row and
#' annotation_col track colors manually. It is  possible to define the colors
#' for only some of the features. Check examples for  details.
#' @param annotation_legend boolean value showing if the legend for annotation
#' tracks should be drawn.
#' @param annotation_names_row boolean value showing if the names for row annotation
#' tracks should be drawn.
#' @param annotation_names_col boolean value showing if the names for column annotation
#' tracks should be drawn.
#' @param drop_levels logical to determine if unused levels are also shown in
#' the legend
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param main the title of the plot
#' @param fontsize base fontsize for the plot
#' @param fontsize_row fontsize for rownames (Default: fontsize)
#' @param fontsize_col fontsize for colnames (Default: fontsize)
#' @param display_numbers logical determining if the numeric values are also printed to
#' the cells. If this is a matrix (with same dimensions as original matrix), the contents
#' of the matrix are shown instead of original values.
#' @param number_format format strings (C printf style) of the numbers shown in cells.
#' For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential
#' notation (see more in \code{\link{sprintf}}).
#' @param number_color color of the text
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param gaps_row vector of row indices that show shere to put gaps into
#'  heatmap. Used only if the rows are not clustered. See \code{cutree_row}
#'  to see how to introduce gaps to clustered rows.
#' @param gaps_col similar to gaps_row, but for columns.
#' @param labels_row custom labels for rows that are used instead of rownames.
#' @param labels_col similar to labels_row, but for columns.
#' @param filename file path where to save the picture. Filetype is decided by
#' the extension in the path. Currently following formats are supported: png, pdf, tiff,
#'  bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is
#' calculated so that the plot would fit there, unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param height manual option for determining the output file height in inches.
#' @param silent do not draw the plot (useful when using the gtable output)
#' @param row_label row cluster labels for semi-clustering
#' @param col_label column cluster labels for semi-clustering
#' @param \dots graphical parameters for the text used in plot. Parameters passed to
#' \code{\link{grid.text}}, see \code{\link{gpar}}.
#'
#' @return
#' Invisibly a list of components
#' \itemize{
#'     \item \code{tree_row} the clustering of rows as \code{\link{hclust}} object
#'     \item \code{tree_col} the clustering of columns as \code{\link{hclust}} object
#'     \item \code{kmeans} the kmeans clustering of rows if parameter \code{kmeans_k} was
#' specified
#' }
#'
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' #@examples
#' # Create test matrix
#' test = matrix(rnorm(200), 20, 10)
#' test[seq(10), seq(1, 10, 2)] = test[seq(10), seq(1, 10, 2)] + 3
#' test[seq(11, 20), seq(2, 10, 2)] = test[seq(11, 20), seq(2, 10, 2)] + 2
#' test[seq(15, 20), seq(2, 10, 2)] = test[seq(15, 20), seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#'
#' # Draw heatmaps
#' pheatmap(test)
#' pheatmap(test, kmeansK = 2)
#' pheatmap(test, scale = "row", clusteringDistanceRows = "correlation")
#' pheatmap(test, color = colorRampPalette(c("navy",
#'     "white", "firebrick3"))(50))
#' pheatmap(test, cluster_row = FALSE)
#' pheatmap(test, legend = FALSE)
#'
#' # Show text within cells
<<<<<<< HEAD
#' pheatmap(test, displayNumbers = TRUE)
#' pheatmap(test, displayNumbers = TRUE, numberFormat = "\%.1e")
#' pheatmap(test, displayNumbers = matrix(ifelse(test > 5,
#'     "*", ""), nrow(test)))
#' pheatmap(test, cluster_row = FALSE,
#'     legendBreaks = seq(-1, 4), legendLabels = c("0",
#'     "1e-4", "1e-3", "1e-2", "1e-1", "1"))
#'
#' # Fix cell sizes and save to file with correct size
#' pheatmap(test, cellWidth = 15, cellHeight = 12, main = "Example heatmap")
#' pheatmap(test, cellWidth = 15, cellHeight = 12, fontSize = 8,
#'     fileName = "test.pdf")
#'
#' # Generate annotations for rows and columns
#' annotationCol = data.frame(CellType = factor(rep(c("CT1", "CT2"), 5)),
#'     Time = seq(5))
#' rownames(annotationCol) = paste("Test", seq(10), sep = "")
#'
#' annotationRow = data.frame(GeneClass = factor(rep(c("Path1",
#'    "Path2",
#'    "Path3"),
#'    c(10, 4, 6))))
#' rownames(annotationRow) = paste("Gene", seq(20), sep = "")
#'
#' # Display row and color annotations
#' pheatmap(test, annotationCol = annotationCol)
#' pheatmap(test, annotationCol = annotationCol, annotationLegend = FALSE)
#' pheatmap(test, annotationCol = annotationCol, annotationRow = annotationRow)
=======
#' pheatmap(test, display_numbers = TRUE)
#' pheatmap(test, display_numbers = TRUE, number_format = "\%.1e")
#' pheatmap(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))
#' pheatmap(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
#' "1e-4", "1e-3", "1e-2", "1e-1", "1"))
#'
#' # Fix cell sizes and save to file with correct size
#' pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
#' pheatmap(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")
#'
#' # Generate annotations for rows and columns
#' annotation_col = data.frame(
#'                     CellType = factor(rep(c("CT1", "CT2"), 5)),
#'                     Time = 1:5
#'                 )
#' rownames(annotation_col) = paste("Test", 1:10, sep = "")
#'
#' annotation_row = data.frame(
#'                     GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
#'                 )
#' rownames(annotation_row) = paste("Gene", 1:20, sep = "")
#'
#' # Display row and color annotations
#' pheatmap(test, annotation_col = annotation_col)
#' pheatmap(test, annotation_col = annotation_col, annotation_legend = FALSE)
#' pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row)
#'
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
#'
#' # Specify colors
#' ann_colors = list(Time = c("white", "firebrick"),
#'     CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
<<<<<<< HEAD
#'     GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))
#'
#' pheatmap(test, annotationCol = annotationCol, annotationColors = ann_colors,
#'     main = "Title")
#' pheatmap(test, annotationCol = annotationCol, annotationRow = annotationRow,
#'     annotationColors = ann_colors)
#' pheatmap(test, annotationCol = annotationCol,
#'     annotationColors = ann_colors[2])
#'
#' # Gaps in heatmaps
#' pheatmap(test, annotationCol = annotationCol, clusterRows = FALSE,
#'     gapsRow = c(10, 14))
#' pheatmap(test, annotationCol = annotationCol, clusterRows = FALSE,
#'     gapsRow = c(10, 14), cutreeCol = 2)
#'
#' # Show custom strings as row/col names
#' labelsRow = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
#' "", "", "Il10", "Il15", "Il1b")
#'
#' pheatmap(test, annotationCol = annotationCol, labelsRow = labelsRow)
=======
#'     GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
#' )
#'
#' pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors, main = "Title")
#' pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row,
#'          annotation_colors = ann_colors)
#' pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors[2])
#'
#' # Gaps in heatmaps
#' pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14))
#' pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14),
#'          cutree_col = 2)
#'
#' # Show custom strings as row/col names
#' labels_row = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
#' "", "", "Il10", "Il15", "Il1b")
#'
#' pheatmap(test, annotation_col = annotation_col, labels_row = labels_row)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
#'
#' # Specifying clustering from distance matrix
#' drows = stats::dist(test, method = "minkowski")
#' dcols = stats::dist(t(test), method = "minkowski")
<<<<<<< HEAD
#' pheatmap(test,
#'     clusteringDistanceRows = drows,
#'     clusteringDistanceCols = dcols)
=======
#' pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
#'
#' # Modify ordering of the clusters using clustering callback option
#' callback = function(hc, mat){
#'     sv = svd(t(mat))$v[, 1]
#'     dend = reorder(as.dendrogram(hc), wts = sv)
#'     as.hclust(dend)
#' }
#'
<<<<<<< HEAD
#' pheatmap(test, clusteringCallback = callback)
=======
#' pheatmap(test, clustering_callback = callback)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
#'
#' dontrun{
#' # Same using dendsort package
#' library(dendsort)
#'
#' callback = function(hc, ...){dendsort(hc)}
#' pheatmap(test, clusteringCallback = callback)
#' }
#'
#' @export
<<<<<<< HEAD
semiPheatmap <- function(mat,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    kmeansK = NA,
    breaks = NA,
    borderColor = "grey60",
    cellWidth = NA,
    cellHeight = NA,
    scale = "none",
    clusterRows = TRUE,
    clusterCols = TRUE,
    clusteringDistanceRows = "euclidean",
    clusteringDistanceCols = "euclidean",
    clusteringMethod = "complete",
    clusteringCallback = .identity2,
    cutreeRows = NA,
    cutreeCols = NA,
    treeHeightRow = ifelse(clusterRows, 50, 0),
    treeHeightCol = ifelse(clusterCols, 50, 0),
    legend = TRUE,
    legendBreaks = NA,
    legendLabels = NA,
    annotationRow = NA,
    annotationCol = NA,
    annotation = NA,
    annotationColors = NA,
    annotationLegend = TRUE,
    annotationNamesRow = TRUE,
    annotationNamesCol = TRUE,
    dropLevels = TRUE,
    showRownames = TRUE,
    showColnames = TRUE,
    main = NA,
    fontSize = 10,
    fontSizeRow = fontSize,
    fontSizeCol = fontSize,
    displayNumbers = FALSE,
    numberFormat = "%.2f",
    numberColor = "grey30",
    fontSizeNumber = 0.8 * fontSize,
    gapsRow = NULL,
    gapsCol = NULL,
    labelsRow = NULL,
    labelsCol = NULL,
    fileName = NA,
    width = NA,
    height = NA,
    silent = FALSE,
    rowLabel,
    colLabel,
    ...) {
=======
semi_pheatmap = function(mat,
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
	kmeans_k = NA,
	breaks = NA,
	border_color = "grey60",
	cellwidth = NA,
	cellheight = NA,
	scale = "none",
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	clustering_distance_rows = "euclidean",
	clustering_distance_cols = "euclidean",
	clustering_method = "complete",
	clustering_callback = identity2,
	cutree_rows = NA,
	cutree_cols = NA,
	treeheight_row = ifelse(cluster_rows, 50, 0),
	treeheight_col = ifelse(cluster_cols, 50, 0),
	legend = TRUE, legend_breaks = NA, legend_labels = NA,
	annotation_row = NA, annotation_col = NA,
	annotation = NA, annotation_colors = NA,
	annotation_legend = TRUE, annotation_names_row = TRUE, annotation_names_col = TRUE,
	drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE,
	main = NA, fontsize = 10,
	fontsize_row = fontsize, fontsize_col = fontsize,
	display_numbers = FALSE,
	number_format = "%.2f",
	number_color = "grey30",
	fontsize_number = 0.8 * fontsize,
	gaps_row = NULL, gaps_col = NULL, labels_row = NULL, labels_col = NULL,
	filename = NA, width = NA, height = NA, silent = FALSE,
    row_label, col_label,
    ...){
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    # Set labels
    if (is.null(labelsRow) & !is.null(rownames(mat))) {
        labelsRow <- rownames(mat)
    }
    if (is.null(labelsRow) & is.null(rownames(mat))) {
        labelsRow <- 1:nrow(mat)
        rownames(mat) <- 1:nrow(mat)
    }

    if (is.null(labelsCol) & !is.null(colnames(mat))) {
        labelsCol <- colnames(mat)
    }
    if (is.null(labelsCol) & is.null(colnames(mat))) {
        labelsCol <- 1:ncol(mat)
        colnames(mat) <- 1:ncol(mat)
    }


<<<<<<< HEAD
    if (.is.na2(breaks)) {
        breaks <- .generateBreaks(mat, length(color), center = TRUE)
=======
    if(is.na2(breaks)){
      breaks = generate_breaks(mat, length(color), center = TRUE)
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
    }


    # Kmeans
    if (!is.na(kmeansK)) {
        # Cluster data
<<<<<<< HEAD
        km <- stats::kmeans(mat, kmeansK, iter.max = 100)
        mat <- km$centers
=======
        km = stats::kmeans(mat, kmeans_k, iter.max = 100)
        mat = km$centers
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

        # Compose rownames
        t <- table(km$cluster)
        labelsRow <- sprintf("Cluster: %s Size: %d", names(t), t)
    } else {
        km <- NA
    }

    # Format numbers to be displayed in cells
<<<<<<< HEAD
    if (is.matrix(displayNumbers) | is.data.frame(displayNumbers)) {
        if (nrow(displayNumbers) != nrow(mat) |
                ncol(displayNumbers) != ncol(mat)) {
            stop("If displayNumbers provided as matrix,
                its dimensions have to match with mat")
=======
    if(is.matrix(display_numbers) | is.data.frame(display_numbers)){
        if(nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)){
            stop("If display_numbers provided as matrix, its dimensions have to match with mat")
        }

        display_numbers = as.matrix(display_numbers)
        fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
        fmat_draw = TRUE

    }
    else{
        if(display_numbers){
            fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
            fmat_draw = TRUE
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
        }

        displayNumbers <- as.matrix(displayNumbers)
        fmat <- matrix(as.character(displayNumbers),
                    nrow = nrow(displayNumbers),
                    ncol = ncol(displayNumbers))
        fmatDraw <- TRUE
    } else {
        if (displayNumbers) {
            fmat <- matrix(sprintf(numberFormat, mat),
                    nrow = nrow(mat),
                    ncol = ncol(mat))
            fmatDraw <- TRUE
        } else {
            fmat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
            fmatDraw <- FALSE
        }
    }

    # Do clustering for rows
<<<<<<< HEAD
    if (clusterRows == TRUE) {
        if (is.null(rowLabel)) {
            rowLabel <- rep(1, nrow(mat))
        } else {
            o <- order(rowLabel)
            mat <- mat[o, , drop = FALSE]
            fmat <- fmat[o, , drop = FALSE]
            rowLabel <- rowLabel[o]
            if (!is.null(annotationRow)) {
                annotationRow <- annotationRow[o, , drop = FALSE]
            }
        }

        treeRow <- .clusterMat(mat,
            rowLabel,
            distance = clusteringDistanceRows,
            method = clusteringMethod)
        treeRow <- clusteringCallback(treeRow, mat)

        mat <- mat[treeRow$order, , drop = FALSE]
        fmat <- fmat[treeRow$order, , drop = FALSE]
        labelsRow <- labelsRow[treeRow$order]
        if (!is.na(cutreeRows)) {
            gapsRow <- .findGaps(treeRow, cutreeRows)
        } else {
            gapsRow <- NULL
        }
    } else {
        treeRow <- NA
        treeHeightRow <- 0
=======
    if(cluster_rows == TRUE) {
      if (is.null(row_label)) {
        row_label = rep(1, nrow(mat))
      } else {
        o = order(row_label)
        mat = mat[o,,drop=FALSE]
        fmat = fmat[o, , drop = FALSE]
        row_label = row_label[o]
        if(!is.null(annotation_row)) {
          annotation_row = annotation_row[o,,drop=FALSE]
        }
      }

   	  tree_row = cluster_mat(mat, row_label , distance = clustering_distance_rows, method = clustering_method)
      tree_row = clustering_callback(tree_row, mat)

      mat = mat[tree_row$order, , drop = FALSE]
      fmat = fmat[tree_row$order, , drop = FALSE]
      labels_row = labels_row[tree_row$order]
      if(!is.na(cutree_rows)){
        gaps_row = find_gaps(tree_row, cutree_rows)
      }
      else {
        gaps_row = NULL
      }
    } else {
      tree_row = NA
      treeheight_row = 0
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
    }


    ## Do clustering for columns
<<<<<<< HEAD
    if (clusterCols == TRUE) {
        if (is.null(colLabel)) {
            colLabel <- rep(1, ncol(mat))
        } else {
            o <- order(colLabel)
            mat <- mat[, o, drop = FALSE]
            fmat <- fmat[, o, drop = FALSE]
            colLabel <- colLabel[o]
            if (!is.null(annotationCol)) {
                annotationCol <- annotationCol[o, , drop = FALSE]
            }
        }

        treeCol <- .clusterMat(t(mat),
            colLabel,
            distance = clusteringDistanceCols,
            method = clusteringMethod)
        treeCol <- clusteringCallback(treeCol, t(mat))

        mat <- mat[, treeCol$order, drop = FALSE]
        fmat <- fmat[, treeCol$order, drop = FALSE]
        labelsCol <- labelsCol[treeCol$order]

        if (!is.na(cutreeCols)) {
            gapsCol <- .findGaps(treeCol, cutreeCols)
        } else {
            gapsCol <- NULL
        }
    } else {
        treeCol <- NA
        treeHeightCol <- 0
    }

    attr(fmat, "draw") <- fmatDraw
=======
    if(cluster_cols == TRUE) {
      if(is.null(col_label)) {
        col_label = rep(1, ncol(mat))
      } else {
        o = order(col_label)
        mat = mat[,o,drop=FALSE]
        fmat = fmat[,o,drop = FALSE]
        col_label = col_label[o]
        if(!is.null(annotation_col)) {
          annotation_col = annotation_col[o,,drop=FALSE]
        }
      }


      tree_col = cluster_mat(t(mat), col_label , distance = clustering_distance_cols, method = clustering_method)
      tree_col = clustering_callback(tree_col, t(mat))

      mat = mat[, tree_col$order, drop = FALSE]
      fmat = fmat[, tree_col$order, drop = FALSE]
      labels_col = labels_col[tree_col$order]

      if(!is.na(cutree_cols)){
        gaps_col = find_gaps(tree_col, cutree_cols)
      }
      else {
        gaps_col = NULL
      }
    } else {
      tree_col = NA
      treeheight_col = 0
    }

    attr(fmat, "draw") = fmat_draw
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    # Colors and scales
    if (!.is.na2(legendBreaks) & !.is.na2(legendLabels)) {
        if (length(legendBreaks) != length(legendLabels)) {
            stop("Lengths of legendBreaks and legendLabels must be the same")
        }
    }


<<<<<<< HEAD
    if (.is.na2(breaks)) {
        breaks <- .generateBreaks(as.vector(mat), length(color))
=======
    if(is.na2(breaks)){
        breaks = generate_breaks(as.vector(mat), length(color))
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
    }
    if (legend & .is.na2(legendBreaks)) {
        legend <- grid.pretty(range(as.vector(breaks)))
        names(legend) <- legend
    }
<<<<<<< HEAD
    else if (legend & !.is.na2(legendBreaks)) {
        legend <- legendBreaks[legendBreaks >= min(breaks) &
                legendBreaks <= max(breaks)]

        if (!.is.na2(legendLabels)) {
            legendLabels <- legendLabels[legendBreaks >= min(breaks) &
                    legendBreaks <= max(breaks)]
            names(legend) <- legendLabels
        } else {
            names(legend) <- legend
=======
    else if(legend & !is.na2(legend_breaks)){
        legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]

        if(!is.na2(legend_labels)){
            legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
            names(legend) = legend_labels
        }
        else{
            names(legend) = legend
>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd
        }
    } else {
        legend <- NA
    }
<<<<<<< HEAD
    mat <- .scaleColours(mat, col = color, breaks = breaks)

    annotation <- c(annotationRow, annotationCol)
    annotation <- annotation[unlist(lapply(annotation,
        function(x) !.is.na2(x)))]
    if (length(annotation) != 0) {
        annotationColors <- .generateAnnotationColours(annotation,
            annotationColors,
            drop = dropLevels)
    } else {
        annotationColors <- NA
    }

    labelsRow <- rownames(mat)
    labelsCol <- colnames(mat)

    if (!showRownames) {
        labelsRow <- NULL
    }

    if (!showColnames) {
        labelsCol <- NULL
    }

    # Draw heatmap
    gt <- .heatmapMotor(mat,
        borderColor = borderColor,
        cellWidth = cellWidth,
        cellHeight = cellHeight,
        treeHeightCol = treeHeightCol,
        treeHeightRow = treeHeightRow,
        treeCol = treeCol,
        treeRow = treeRow,
        fileName = fileName,
        width = width,
        height = height,
        breaks = breaks,
        color = color,
        legend = legend,
        annotationRow = annotationRow,
        annotationCol = annotationCol,
        annotationColors = annotationColors,
        annotationLegend = annotationLegend,
        annotationNamesRow = annotationNamesRow,
        annotationNamesCol = annotationNamesCol,
        main = main,
        fontSize = fontSize,
        fontSizeRow = fontSizeRow,
        fontSizeCol = fontSizeCol,
        fmat = fmat,
        fontSizeNumber = fontSizeNumber,
        numberColor = numberColor,
        gapsRow = gapsRow,
        gapsCol = gapsCol,
        labelsRow = labelsRow,
        labelsCol = labelsCol,
        ...)

    if (is.na(fileName) & !silent) {
        grid.newpage()
        grid.draw(gt)
    }
=======
    else {
        legend = NA
    }
    mat = scale_colours(mat, col = color, breaks = breaks)

    annotation = c(annotation_row, annotation_col)
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]
    if(length(annotation) != 0){
        annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
    }
    else{
        annotation_colors = NA
    }

    labels_row = rownames(mat)
    labels_col = colnames(mat)

    if(!show_rownames){
        labels_row = NULL
    }

    if(!show_colnames){
        labels_col = NULL
    }

    # Draw heatmap
    gt = heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, number_color = number_color, gaps_row = gaps_row, gaps_col = gaps_col, labels_row = labels_row, labels_col = labels_col, ...)

    if(is.na(filename) & !silent){
        grid.newpage()
        grid.draw(gt)
    }

    invisible(list(tree_row = tree_row, tree_col = tree_col, gtable = gt))
}

>>>>>>> c9e56fdfe5198b640f49fe8f6c2b203e37a2d3fd

    invisible(list(treeRow = treeRow,
        treeCol = treeCol,
        gtable = gt))
}
