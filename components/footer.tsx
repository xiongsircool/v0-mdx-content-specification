import Link from "next/link"
import { Dna, Github, Mail, MapPin } from "lucide-react"

export function Footer() {
  return (
    <footer className="border-t bg-muted/50">
      <div className="container mx-auto px-4 py-12">
        <div className="grid grid-cols-1 md:grid-cols-4 gap-8">
          {/* Logo and Description */}
          <div className="md:col-span-2">
            <div className="flex items-center gap-3 mb-4">
              <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-primary text-primary-foreground">
                <Dna className="h-5 w-5" />
              </div>
              <div>
                <h3 className="text-lg font-bold">SBC</h3>
                <p className="text-sm text-muted-foreground">上海生物信息中心学生组织</p>
              </div>
            </div>
            <p className="text-sm text-muted-foreground text-pretty max-w-md">
              致力于推广生物信息学知识，为学生提供学习资源、技术分享和学术交流平台。
              我们相信通过协作和知识共享，能够推动生物信息学领域的发展。
            </p>
          </div>

          {/* Quick Links */}
          <div>
            <h4 className="text-sm font-semibold mb-4">快速导航</h4>
            <ul className="space-y-2 text-sm">
              <li>
                <Link href="/announcements" className="text-muted-foreground hover:text-foreground transition-colors">
                  最新公告
                </Link>
              </li>
              <li>
                <Link href="/posts" className="text-muted-foreground hover:text-foreground transition-colors">
                  技术推文
                </Link>
              </li>
              <li>
                <Link href="/learn" className="text-muted-foreground hover:text-foreground transition-colors">
                  学习资源
                </Link>
              </li>
              <li>
                <Link href="/about" className="text-muted-foreground hover:text-foreground transition-colors">
                  关于我们
                </Link>
              </li>
            </ul>
          </div>

          {/* Contact */}
          <div>
            <h4 className="text-sm font-semibold mb-4">联系我们</h4>
            <ul className="space-y-2 text-sm">
              <li className="flex items-center gap-2 text-muted-foreground">
                <Mail className="h-4 w-4" />
                <span>contact@sbc.org.cn</span>
              </li>
              <li className="flex items-center gap-2 text-muted-foreground">
                <Github className="h-4 w-4" />
                <a href="#" className="hover:text-foreground transition-colors">
                  GitHub
                </a>
              </li>
              <li className="flex items-center gap-2 text-muted-foreground">
                <MapPin className="h-4 w-4" />
                <span>上海市</span>
              </li>
            </ul>
          </div>
        </div>

        <div className="border-t mt-8 pt-8 flex flex-col sm:flex-row justify-between items-center gap-4">
          <p className="text-xs text-muted-foreground">© 2025 上海生物信息中心学生组织. 保留所有权利.</p>
          <div className="flex gap-4 text-xs text-muted-foreground">
            <Link href="/privacy" className="hover:text-foreground transition-colors">
              隐私政策
            </Link>
            <Link href="/terms" className="hover:text-foreground transition-colors">
              使用条款
            </Link>
          </div>
        </div>
      </div>
    </footer>
  )
}
