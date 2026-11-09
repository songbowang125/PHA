import os
import multiprocessing


def download_one(output_file, line):
    os.system("wget -c -O {} {}".format(output_file, line))


if __name__ == '__main__':


    input_file = "download_files_part1"

    pool = multiprocessing.Pool(10)

    with open(input_file) as fin:
        for line in fin:

            line = line.strip()

            file_name = os.path.basename(line)
            sample_name = line.split("/")[-2]

            output_file = "{}.{}".format(sample_name, file_name)

            pool.apply_async(download_one, (output_file, line))


    pool.close()
    pool.join()