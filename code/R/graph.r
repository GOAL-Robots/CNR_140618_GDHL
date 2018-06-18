#!/usr/local/bin/Rexec
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# usage: graph.r <name>    name of the dataset
#                <nvars>   number of parameters used
#                <index>   index of the parameters' permutation
#                <seed>    seed of the genetic algorithm
#                [scale]   scaling of original data [default to 1]
#                [trans]   translation of original data [default to 0]
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Init arguments list ----
args <- commandArgs()
args_start = grep("--args", args) + 1

if (length(args_start) > 0 && args_start < length(args)) {
    args <- args[args_start:length(args)]
} else {
    args = c()
}

require(data.table)

# Utilities ----

idxparms <- function(nvars, index) {
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #  Each index is a specific permutation of
    #  the current number of parameters used over
    #  the 8 possible parameters.
    #
    #  example
    #  nvars = 5; index = 3
    #  permutation = 0 0 1 1 0 1 1 1
    #
    #  that is, simulations with nvars(5
    #  and index(3) will use the 3rd, 4th
    #  6th, 7th and 8th kernel of the
    #  equation, and only the corresponding
    #  parameters will be searched for.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    require(gtools)



    p = permutations(8, 8)

    # put to zero parameters > nvars
    p = p * (p > (8 - nvars)) - (8 - nvars)

    p = p * (p > 0)

    # put to one all other parameters
    # we need only the position in
    # the sequence
    p = unique(1 * (p > 0))

    # get the sequence at
    # position "index"
    p = p[index,]

    return(p)
}

g_table <- function(nvars, params, FVU, file_name) {
    tabl = c()

    tabl = rbind(tabl, c('Number of components:', nvars))
    tabl = rbind(tabl, c('FVU:', FVU))
    tabl = rbind(tabl, c('Onset amplification:', params[3]))
    tabl = rbind(tabl, c('pre-synaptic tau:', params[1]))
    tabl = rbind(tabl, c('post-synaptic tau:', params[2]))

    require(xtable)
    print.xtable(
        xtable(tabl, align = "rrl"),
        type = "latex",
        file = paste(file_name, '-params_table.tex', sep = ''),
        include.rownames = FALSE,
        include.colnames = FALSE,
        hline.after = FALSE
    )
}

g_kernels <- function(pars, file_name) {
    postscript(
        paste(file_name, '-used_kernels.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        pointsize=6,
        width = 2.5,
        height = 1.67
    )

    par(mar = c(5, 4, 2, 2) + .1)
    pidx = (1:8)[pars != 0]

    labels = c(
        expression(sigma[pp]),
        expression(sigma[np]),
        expression(sigma[pn]),
        expression(sigma[nn]),
        expression(eta[ps]),
        expression(eta[sp]),
        expression(eta[ns]),
        expression(eta[sn])
    )

    plot(
        pidx,
        pars[pidx],
        cex = 1,
        cex.lab = 2,
        xlim = c(.5, 8.5),
        ylim = c(-1.3, 1.3),
        axes = FALSE,
        xlab = "Component type",
        ylab = ''
    )

    pp = pars[pidx]
    pp = as.matrix(pp)
    arrows(pidx,
           pp * 0,
           pidx,
           pp,
           lwd = 1,
           length = 0.08)

    points(1:8, pars * 0, cex = .05)


    for (x in 1:length(pidx)) {
        text(pidx[x],
             pars[pidx[x]]
             + .3 * (pars[pidx[x]] / abs(pars[pidx[x]])),
             format(pars[pidx[x]],
                    digits = 2),
             cex = 1.2)
    }

    axis(
        1,
        at = 1:8,
        labels = labels,
        lwd = 0,
        cex.axis = 1.5
    )

    axis(2,
         at = seq(-1, 1, 1),
         lwd = 1,
         cex.axis = 1.5)

    dev.off()

}

g_fitting <- function(data, regression, file_name) {
    postscript(
        paste(file_name, '-data_fitting.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        width = 6,
        height = 6,

    )


    data = data[order(data[, 1]),]
    data_plus_regression = c(data[, 2], regression[, 2])

    dr = c(max(data[, 1]) - min(data[, 1]),
           max(data_plus_regression) - min(data_plus_regression))
    drg = c(dr[1] * (.1 / 2), dr[2] * (.1 / 2))

    ylim = c(min(data_plus_regression) - drg[2] ,
             max(data_plus_regression) + drg[2])

    xlim = c(min(data[, 1]) - drg[1] ,  max(data[, 1]) + drg[1])

    print(xlim)
    plot(
        axes = FALSE,
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        data,
        ylim = ylim,
        xlim = xlim
    )

    arrows(
        c(min(data[, 1]) - drg[1]  ,  0),
        c(0                       ,  min(data_plus_regression) - drg[2]),
        c(max(data[, 1]) + drg[1]  ,  0),
        c(0                       ,  max(data_plus_regression) + drg[2]),
        length = .08,
        lwd = 3
    )

    lines(regression, lwd = 6)
    points(data, pch = 19, col = "#ffffff")

    points(data)


    lab = "Delta*t"
    lab = eval(parse(text = paste('expression(', lab, ')')))

    text(xlim[2],-.1 * (ylim[2] - ylim[1]), lab, cex = 2)
    lab = "-Delta*t"
    lab = eval(parse(text = paste('expression(', lab, ')')))

    text(xlim[1],-.1 * (ylim[2] - ylim[1]), lab, cex = 2)

    axis(
        1,
        lwd = 0,
        lwd.ticks = .1,
        at = seq(-100, 100, 50),
        cex.axis = 1.5
    )

    dev.off()
}

g_clear <- function(data, file_name, data_name) {
    postscript(
        paste(file_name, '-clear.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        pointsize = 12,
        width = 6,
        height = 6
    )

    load("data_bezier.dat")

    labels = list()
    labels[["BiPoo2001"]]                  = list(
        "x" = "Spike timing (msec)",
        "y" = "Synaptic change (%)",
        "atx" = seq(120,-100,-40),
        "aty" = seq(-60, 100, 20),
        "laby" = seq(-60, 100, 20)
    )
    labels[["WoodinGangulyPoo2003"]]       = list(
        "x" = "Spike timing (msec)",
        "y" = "Change in GPSC amplitude (%)",
        "atx" = seq(-80, 80, 40),
        "aty" = seq(-50, 100, 50),
        "laby" = seq(-50, 100, 50)
    )
    labels[["CassenaerLaurent2007"]]       = list(
        "x" = "Spike timing (ms)",
        "y" = "Synaptic change (%)",
        "atx" = seq(-40, 40, 20),
        "aty" = seq(-60, 100, 20),
        "laby" = seq(-60, 100, 20)
    )
    labels[["FroemkeDan2002"]]             = list(
        "x" = "Spike timing (msec)",
        "y" = "Synaptic change (%)",
        "atx" = seq(-40, 40, 20),
        "aty" = seq(-60, 100, 20),
        "laby" = seq(-60, 100, 20)
    )
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    labels[["ZhouAckerWhite2005"]]         = list(
        "x" = "Spike timing (ms)",
        "y" = "% Change in EPSPs",
        "atx" = seq(-40, 40, 20),
        "aty" = seq(-80, 80, 40),
        "laby" = seq(-80, 80, 40)
    )
    labels[["ZhouAckerWhite2005-0"]]       = list(
        "x" = "Spike timing (ms)",
        "y" = "% Change in EPSPs",
        "atx" = seq(-40, 40, 20),
        "aty" = seq(-60, 100, 20),
        "laby" = seq(-60, 100, 20)
    )
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    labels[["HaasNowotnyAbarbanel2006"]]   = list(
        "x" = expression(paste(Delta, "t [ms]")),
        "y" = expression(paste(Delta[ISPS], "[norm]")),
        "atx" = seq(-40, 40, 20),
        "aty" = seq(-50, 100, 50),
        "laby" = seq(0.5, 2, 0.5)
    )
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    labels[["WittenbergWang20062"]]        = list(
        "x" = "Delay from EPSP to 2nd AP (mc)",
        "y" = "Normalized synaptic strength",
        "atx" = seq(-80, 80, 40),
        "aty" = seq(-50, 100, 50),
        "laby" = seq(0.5, 2, 0.5)
    )
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    labels[["ZhangTaoPoo1998"]]            = list(
        "x" = "Spike timing (msec)",
        "y" = "Synaptic change (%)",
        "atx" = seq(-40, 40, 20),
        "aty" = seq(-60, 100, 20),
        "laby" = seq(-60, 100, 20)
    )
    par(mar = c(5, 6, 2, 2) + .1)

    data = data[order(data[, 1]),]

    regression = data
    regression[, 2] = 0
    data_plus_regression = c(data[, 2], regression[, 2])

    dr = c(max(data[, 1]) - min(data[, 1]),
           max(data_plus_regression) - min(data_plus_regression))
    drg = c(dr[1] * (.1 / 2), dr[2] * (.1 / 2))

    ylim = c(min(min(data_plus_regression) - drg[2], labels[[data_name]][["aty"]])  ,
             max(max(data_plus_regression) + drg[2], labels[[data_name]][["aty"]]))

    xlim = c(min(data[, 1]) - drg[1] ,  max(data[, 1]) + drg[1])

    print(xlim)
    plot(
        axes = FALSE,
        frame.plot = TRUE,
        xlab = labels[[data_name]][["x"]],
        ylab = labels[[data_name]][["y"]],
        cex.lab = 2,
        data,
        ylim = ylim,
        xlim = xlim
    )

    if (data_name %in% data_bezier$name) {
        curr_bezier = data_bezier[data_bezier$name == data_name]
        bez = curr_bezier[idx == 1]
        lines(bez$x, bez$y, lwd = 3, col = "#ff0000")
        bez = curr_bezier[idx == 2]
        lines(bez$x, bez$y, lwd = 3, col = "#ff0000")
    }

    arrows(
        c(min(data[, 1]) - drg[1]  ,  0),
        c(0                       ,  min(data_plus_regression) - drg[2]),
        c(max(data[, 1]) + drg[1]  ,  0),
        c(0                       ,  max(data_plus_regression) + drg[2]),
        length = .08,
        lwd = 3
    )

    points(data, pch = 19, col = "#ffffff")

    points(data)

    axis(
        1,
        lwd = .0,
        lwd.ticks = .5,
        at = labels[[data_name]][["atx"]],
        cex.axis = 1.5
    )

    axis(
        2,
        lwd = .0,
        lwd.ticks = .5,
        at = labels[[data_name]][["aty"]],
        labels = labels[[data_name]][["laby"]],
        cex.axis = 1.5
    )
    dev.off()
}

g_bic <- function(bics, current, file_name) {
    postscript(
        paste(file_name, '-number_of_kernels.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        pointsize=6,
        width = 2.5,
        height = 1.67
    )

    par(mar = c(5, 6, 2, 2) + .1)

    bicrange = c(min(bics) - abs(max(bics) - min(bics)) / 6 ,
                 max(bics) + abs(max(bics) - min(bics)) / 6)
    plot(
        bics,
        pch = 15,
        frame.plot = FALSE,
        axes = FALSE,
        ylim = bicrange,
        xlab = "Number of components",
        ylab = "BIC",
        cex = 1.5,
        cex.lab = 2,
        cex.axis = 2
    )

    lines(bics, lwd = 4)

    points(
        which.min(bics),
        min(bics),
        pch = 1,
        cex = 3.5,
        lwd = 1
    )
    points(current,
           bics[current],
           pch = 2,
           cex = 3.5,
           lwd = 1)

    q = as.integer(bics)
    axis(2, at = seq(min(q), max(q), abs(max(q) - min(q)) / 2) , lwd = 1)
    axis(1, at = seq(1, 8, 1), lwd = 1)

    deltabic =  abs(min(bics) - bics[current])
    if (deltabic > 0)
    {
        lab = "Delta*BIC"
        lab = eval(parse(text = paste(
            'expression(',
            'Delta*BIC',
            '==',
            format(deltabic),
            ')'
        )))

        text(4, max(bics) + abs(max(bics) - min(bics)) / 12, lab, cex =
                 1.5)
    }
    dev.off()


}

# Functions ----------------------------------------------------------------------------------------

graph <- function(name,
                  nvars,
                  index,
                  seed,
                  data_dir,
                  scale = 1,
                  trans = 0) {
    # radix of the data filenames for this kernel combination
    radix_file = paste(
        'curve-',
        name,
        '-',
        format(nvars),
        '-',
        format(index),
        '-',
        format(seed),
        '_',
        sep = ''
    )


    label_file = paste(radix_file,
                       'data',
                       sep = '')


    # radix of the data filenames for this kernel combination
    # best fitting
    curve_file = paste(data_dir, '/', radix_file, 'data', sep = '')

    # experimental_data
    data_file = paste(data_dir, '/', radix_file, 'orig_data', sep = '')

    # parameters for the best fitting
    params_file = paste(data_dir, '/', radix_file, 'params', sep = '')

    # bics
    bics_file = paste(data_dir, '/', 'bics', sep = '')


    data_all = read.table('data_all')

    colnames(data_all) = c(
        'name',
        'nvars',
        'index',
        'seed',
        'FVU',
        'tau1',
        'tau2',
        'k',
        'gpp',
        'gnp',
        'gpn',
        'gnn',
        'eps',
        'esp',
        'ens',
        'esn'
    )

    FVU = data_all[data_all$name == name &
                       data_all$nvars == nvars &
                       data_all$index == index &
                       data_all$seed == seed, ]$FVU


    curve_ = as.matrix(read.table(curve_file))

    data_ = as.matrix(read.table(data_file))

    data_[, 2] = data_[, 2] * scale + trans
    bics_ = read.table(bics_file)

    bics_ = bics_[bics_[, 1] == name,]

    bics_ = bics_[order(bics_[, 2]),][, 5]

    params_ = as.matrix(read.table(params_file))



    # Parameters table ----
    g_table(nvars, params_, FVU, label_file)

    # Used kernels ----
    g_kernels(params_[4:11], label_file)

    # Fitting ----
    g_fitting(data_, curve_, label_file)
    g_clear(data_, label_file, name)

    # BIC -----
    g_bic(bics_, nvars, label_file)

}

custom_graph <- function(dt = .01,
                         tau1 = 1,
                         tau2 = 1,
                         k = 1,
                         gpp = 0,
                         gnp = 0,
                         gpn = 0,
                         gnn = 0,
                         eps = 0,
                         esp = 0,
                         ens = 0,
                         esn = 0,
                         width = 6,
                         height = 4) {
    print(args)
    if (length(args) != 7) {
        return('wrong number of args')

    }

    name = args[1]

    nvars = as.numeric(args[2])

    index = as.numeric(args[3])

    seed = as.numeric(args[4])

    data_dir = args[5]

    scale = as.numeric(args[6])

    trans = as.numeric(args[7])




    curve_ = as.matrix(read.table(curve_file))

    data_ = as.matrix(read.table(data_file))

    data_[, 2] = data_[, 2] * scale + trans
    scores_ = as.matrix(read.table(scores_file))

    params_ = as.matrix(read.table(params_file))




    # Parameters table ----
    tabl = c()

    tabl = rbind(tabl, 'Data:', name)

    tabl = rbind(tabl, 'Kernels number:', nvars)
    tabl = rbind(tabl, 'Seed:', seed)
    tabl = rbind(tabl, 'kernel amplification:', params_[3])
    tabl = rbind(tabl, 'pre-synaptic tau:', params_[1])
    tabl = rbind(tabl, 'post-synaptic tau:', params_[2])

    require(xtable)
    tabl = print(xtable(tabl))
    write(tabl, paste(curve_file, '-params_table.tex', sep = ''))

    # Used kernels ----
    postscript(
        paste(curve_file, '-used_kernels.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        pointsize = 11,
        width = 2.5,
        height = 1.67,

    )

    par(mar = c(5, 4, 2, 2) + .1)

    pars = params_[4:11]

    pidx = (1:8)[pars != 0]

    labels = c(
        expression(sigma[pp]),
        expression(sigma[np]),
        expression(sigma[pn]),
        expression(sigma[nn]),
        expression(eta[ps]),
        expression(eta[sp]),
        expression(eta[ns]),
        expression(eta[sn])
    )



    plot(
        pidx,
        pars[pidx],
        cex = .5,
        cex.lab = 2,
        xlim = c(.5, 8.5),
        ylim = c(-1.3, 1.3),
        axes = FALSE,
        xlab = "GDHL Component type",
        ylab = ''
    )


    pp = pars[pidx]
    pp = as.matrix(pp)
    arrows(pidx,
           pp * 0,
           pidx,
           pp,
           lwd = 4,
           length = 0.08)


    points(1:8, pars * 0, cex = .05)


    for (x in 1:length(pidx)) {
        text(pidx[x],
             pars[pidx[x]]
             + .3 * (pars[pidx[x]] / abs(pars[pidx[x]])),
             format(pars[pidx[x]],
                    digits = 2),
             cex = 1.3)

    }

    axis(
        1,
        at = 1:8,
        labels = labels,
        lwd = 0,
        cex.axis = 2
    )

    axis(2,
         at = seq(-1, 1, 1),
         lwd = 3,
         cex.axis = 1)


    dev.off()

    # GA history ----
    postscript(
        paste(curve_file, '-fitness_curve.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        width = 6,
        height = 6,

    )


    par(mar = c(5, 4, 2, 1) + .1)

    plot(
        1:length(scores_),
        main = "Fitness curve",
        ylab = "Best scores",
        scores_,
        cex = .1,

        xlab = 'Generations',
    )

    grid()


    dev.off()

    # Fitting ----
    postscript(
        paste(curve_file, '-data_fitting.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        width = 6,
        height = 6,

    )

    data_ = data_[order(data_[, 1]),]

    all_scores = c(data_[, 2], curve_[, 2])

    # for(f in similar_files)  {
    # similar_=as.matrix(read.table(f));
    # all_scores=c(all_scores,similar_[,2]);
    # }

    dr = c(max(data_[, 1]) - min(data_[, 1]),
           max(all_scores) - min(all_scores))
    drg = c(dr[1] * (.1 / 2), dr[2] * (.1 / 2))



    ylim = c(min(all_scores) - drg[2] , max(all_scores) + drg[2])

    xlim = c(min(data_[, 1]) - drg[1] ,  max(data_[, 1]) + drg[1])


    plot(
        axes = FALSE,
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        data_,
        ylim = ylim,
        xlim = xlim
    )



    arrows(
        c(min(data_[, 1]) - drg[1], 0),
        c(0,                           min(all_scores) - drg[2]),
        c(max(data_[, 1]) + drg[1], 0),
        c(0,                           max(all_scores) + drg[2]),
        length = .08,
        lwd = 3
    )



    # for(f in similar_files)  {
    # similar_=as.matrix(read.table(f));
    # lines(similar_,col='#aaaaaa',lwd=2);
    # }
    lines(curve_, lwd = 6)

    points(data_, pch = 19, col = "#ffffff")

    points(data_)

    #grid();

    dev.off()

    postscript(
        paste(curve_file, '-number_of_kernels.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        pointsize = 11,
        width = 2.5,
        height = 1.67,

    )


    par(mar = c(5, 6, 2, 2) + .1)


    allscores = c()
    for (f in similar_scores) {
        sc = as.matrix(read.table(f))
        sc = sc[dim(sc)[1],]
        allscores = c(allscores, sc)

    }

    plot(
        allscores,
        pch = 15,
        ylim = c(0, .8),
        frame.plot = FALSE,
        axes = FALSE,
        xlab = "Number of components",
        ylab = "Best fitting",
        cex = 1.5,
        cex.lab = 2,
        cex.axis = 2
    )

    lines(allscores, lwd = 4)

    points(nvars,
           allscores[nvars],
           pch = 1,
           cex = 3.5,
           lwd = 3)

    axis(2, at = seq(0, 1, .4), lwd = 3)
    axis(1, at = seq(1, 8, 1), lwd = 3)

    dev.off()

    postscript(
        paste(curve_file, '_legend_all.eps', sep = ''),
        onefile = FALSE,
        horizontal = FALSE,
        width = 6,
        height = 6,

    )


    layout(matrix(c(1, 5, 2, 2, 1, 5, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4), 4, 4))

    par(mar = c(1, 4, 0, 2) + .1)
    plot.new()

    text(.5, .7, format(paste('Kernels number:', nvars), width = 12))
    text(.5, .4, format(paste('kernel amplification:', params_[3]), width =
                            12))
    text(.5, .3, format(paste('pre-synaptic tau:', params_[1]), width = 12))
    text(.5, .2, format(paste('post-synaptic tau:', params_[2]), width = 12))

    par(mar = c(5, 4, 2, 2) + .1)

    pars = params_[4:11]

    pidx = (1:8)[pars != 0]

    labels = c(
        expression(sigma[pp]),
        expression(sigma[np]),
        expression(sigma[pn]),
        expression(sigma[nn]),
        expression(eta[ps]),
        expression(eta[sp]),
        expression(eta[ns]),
        expression(eta[sn])
    )



    plot(
        pidx,
        pars[pidx],
        main = "Used kernels",
        cex = .1,
        xlim = c(.5, 8.5),
        ylim = c(-1.3, 1.3),
        axes = FALSE,
        xlab = '',
        ylab = ''
    )


    pp = pars[pidx]
    pp = as.matrix(pp)
    arrows(pidx,
           pp * 0,
           pidx,
           pp,
           length = 0.03)


    points(1:8, pars * 0, cex = .05)


    for (x in 1:length(pidx)) {
        text(pidx[x],
             pars[pidx[x]]
             + .15 * (pars[pidx[x]] / abs(pars[pidx[x]])),
             format(pars[pidx[x]],
                    digits = 2),
             cex = .8)

    }

    axis(1,
         at = 1:8,
         labels = labels,
         lwd = 0)

    axis(2, at = seq(-1, 1, 1))

    par(mar = c(5, 4, 2, 1) + .1)

    plot(
        1:length(scores_),
        main = "Fitness curve",
        ylab = "Best scores",
        scores_,
        cex = .1,
        xlab = 'Generations',
    )

    grid()

    data_ = data_[order(data_[, 1]),]


    all_scores = c(data_[, 2], curve_[, 2])

    for (f in similar_files)  {
        similar_ = as.matrix(read.table(f))

        all_scores = c(all_scores, similar_[, 2])

    }

    ylim = c(min(all_scores), max(all_scores))


    plot(
        data_,
        main = "Data fitting",
        ylim = ylim,
        ylab = "EPSP",
        xlab = "time-window (msec)"
    )



    for (f in similar_files)  {
        similar_ = as.matrix(read.table(f))

        lines(similar_, col = '#dadada')

    }
    points(data_)

    lines(curve_)

    grid()

    par(mar = c(5, 4, 2, 2) + .1)

    allscores = c()
    for (f in similar_scores) {
        sc = as.matrix(read.table(f))
        sc = sc[dim(sc)[1],]
        allscores = c(allscores, sc)

    }

    plot(
        allscores,
        pch = 14,
        cex = .5,
        ylim = c(0, 1),
        main = "Current kernel number",
        xlab = "Number of kernels",
        ylab = "Best score"
    )

    lines(allscores)

    points(nvars, allscores[nvars], pch = 1, cex = 1.5)
    dev.off()
}


# MAIN ----
if (length(args) != 7) {
    return('wrong number of args')

} else {
    print(args)
    graph(
        name = args[1],
        nvars = as.numeric(args[2]),
        index = as.numeric(args[3]),
        seed = as.numeric(args[4]),
        data_dir = args[5],
        scale = as.numeric(args[6]),
        trans = as.numeric(args[7])
    )

}
