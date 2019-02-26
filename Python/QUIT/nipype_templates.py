############################ <my_interface> ############################


class QiXXInputSpec(QUITCommandInputSpec):
    # Inputs
    # If it is only one input file call it in_file
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Input file')

    # A lot of QUIT tools have parameater files as either JSON or dict file
    param_file = File(desc='.json file', position=-1, argstr='--json=%s',
                      xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='Param Dict', position=-1,
                             argstr='', mandatory=True, xor=['param_file'])

    # Output file is sometimes necessary to be specified as well
    out_file = File(exists=False, argstr='%s', mandatory=True,
                    position=-2, desc='Output file', gen_file=True)

    # Options (Collection of arguments commonly used between the tools. For consistency)
    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    mask_file = File(
        desc='Only process voxels within the mask', argstr='--mask=%s')


class QiXXOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output file")


class QiXX(QUITCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.<module> import <myInterface>
    >>> interface = <myInterface>(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = '<myInterface>'
    input_spec = QiXXInputSpec
    output_spec = QiXXOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()

        outputs['output_img'] = os.path.abspath(
            self._add_prefix(self.inputs.out_file))

        return outputs


######################################################################

# Input/Output/Option traits:
bool_trait = traits.Bool(desc='Description', argstr='-X')

int_trait = traits.Int(desc='Description', argstr='-X=%d')

string_trait = traits.String(desc='Description', argstr='-X=%s')

# Use exist=True if input. False if output
file_Trait = File(desc='Description', argstr='-X=%s', exists=True)

enum_trait = traits.Enum("Opt1", "Opt2", "Opt3",
                         desc="Description", argstr="-X=%d")

flot_trait = traits.Float(desc='Description', argstr='-X=%f')

dict_trait = traits.Dict(desc='Description', position=1,
                         argstr='', mandatory=True)

# If we have mutually exclusive arguments such as dict and json input for params we use xor
param_dict = traits.Dict(desc='dictionary trait', position=1,
                         argstr='', mandatory=True, xor=['param_file'])
param_file = File(desc='Parameter .json file', position=1, argstr='--json=%s',
                  xor=['param_dict'], mandatory=True, exists=True)
