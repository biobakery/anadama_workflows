import mimetypes

def cat(infiles_list, guess_from=None):
    if not guess_from:
        guess_from = infiles_list[0]

    maj_file_type, min_file_type = mimetypes.guess_type(guess_from)
    if min_file_type == 'gzip':
        return "zcat " + " ".join(infiles_list)
    elif min_file_type == 'bzip2':
        return "bzcat " + " ".join(infiles_list)

    # if it's completely unrecognized, just return cat
    return "cat " + " ".join(infiles_list)

            
