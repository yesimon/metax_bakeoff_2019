
S3_ROOT = config.get('S3_ROOT', '')

ABATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID')
if ABATCH_JOB_ID is not None:
    ARRAY_INDEX = os.environ.get('AWS_BATCH_JOB_ARRAY_INDEX')
    unique_dir = ABATCH_JOB_ID
    if ARRAY_INDEX is not None:
        ARRAY_INDEX = int(ARRAY_INDEX)
        samples = [samples[ARRAY_INDEX]]
        unique_dir = os.path.join(unique_dir, str(ARRAY_INDEX))
        workdir: unique_dir

def is_nonstr_iterable(x):
    '''Tests whether `x` is an Iterable other than a string'''
    return isinstance(x, collections.Iterable) and not isinstance(x,str)


import botocore
import botocore.session
def objectify_remote(file_address, *args, **kwargs):
    if file_address is None:
        raise IOError("%s is None" % file_address)

    # if this is a string, make it a list, otherwise use it as an iterable
    file_list = file_address if is_nonstr_iterable(file_address) else [file_address]

    for index, uri in enumerate(file_list):
        if uri.lower().startswith('s3://'):
            import snakemake.remote.S3

            # if botocore cannot find credentials, try connecting unsigned.
            # This should work for anonymous S3 resources
            # and testing on Travis where no credentials are set.
            # This can be removed if/when Snakemake does the same check itself
            if botocore.session.get_session().get_credentials():
                remote_provider = snakemake.remote.S3.RemoteProvider()
            else:
                remote_provider = snakemake.remote.S3.RemoteProvider(config=botocore.client.Config(signature_version=botocore.UNSIGNED))

            file_list[index] = remote_provider.remote(uri, *args, **kwargs)
        elif uri.lower().startswith('gs://'):
            import snakemake.remote.GS
            remote_provider = snakemake.remote.GS.RemoteProvider()
            file_list[index] = remote_provider.remote(uri, *args, **kwargs)
        elif uri.lower().startswith('sftp://'):
            import snakemake.remote.SFTP
            remote_provider = snakemake.remote.SFTP.RemoteProvider()
            file_list[index] = remote_provider.remote(uri, *args, **kwargs)

    return file_list
