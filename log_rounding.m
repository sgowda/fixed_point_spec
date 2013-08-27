function [] = log_rounding(blkname, loss)
f = fopen('/nas/users/sgowda/rounding_loss', 'a');
fprintf(f, '%s: %g\n', blkname, loss);
fclose(f);
