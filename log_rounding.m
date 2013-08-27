function [] = log_rounding(blkname, orig_bitwidth, new_bitwidth, loss)
f = fopen('/nas/users/sgowda/rounding_loss', 'a');
fprintf(f, '%s: originally %d bits, now %d bits, loss = %g\n', blkname, orig_bitwidth, new_bitwidth, loss);
fclose(f);

