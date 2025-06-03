# Greco Examples

Example projects demonstrating how to use the Greco libraries and circuits.

## Examples

This directory will contain example projects showing how to use the Greco libraries and circuits. Each example will be a standalone Noir project with its own README and documentation.

### Prerequisites

1. Install Noir using noirup:

```bash
curl -L https://raw.githubusercontent.com/noir-lang/noirup/refs/heads/main/install | bash
noirup
```

2. Install Barretenberg (BB) proving backend:

```bash
curl -L https://raw.githubusercontent.com/AztecProtocol/aztec-packages/refs/heads/master/barretenberg/bbup/install | bash
bbup
```

After installation, verify both tools are properly installed:

```bash
nargo --version
bb --version
```

## Contributing

To add a new example:

1. Create a new directory under `examples/`
2. Include a descriptive README
3. Ensure the example is self-contained
4. Add appropriate test cases

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.
