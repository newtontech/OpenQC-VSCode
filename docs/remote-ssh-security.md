# Remote SSH Security

OpenQC remote compute integration uses Paramiko and verifies server host keys by default.
Unknown host keys are rejected unless `allow_unknown_host=True` is passed when creating
`SSHHandler`.

Recommended setup:

1. Add trusted compute hosts to the user's OpenSSH `known_hosts` file with `ssh-keyscan`
   or by connecting once with the OpenSSH client after verifying the fingerprint out of band.
2. Pass `known_hosts_path` when the integration should use a project- or environment-specific
   trust store instead of only the user's default SSH trust store.
3. Use key-based authentication or environment-managed credentials. Do not hardcode passwords
   in project files.
4. Prefer `SSHHandler.execute_args([...])` for remote commands. It shell-quotes arguments before
   sending the command string to the SSH server.

`allow_unknown_host=True` should only be used for controlled bootstrapping flows where the host
fingerprint has been verified through another trusted channel.
